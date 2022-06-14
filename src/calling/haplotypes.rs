use crate::model::{AlleleFreq, Data, Likelihood, Marginal, Posterior, Prior};
use anyhow::Result;
use bio::stats::{bayesian::model::Model, probs::LogProb, PHREDProb, Prob};
use bv::BitVec;
use derefable::Derefable;
use derive_builder::Builder;
use derive_deref::DerefMut;
use derive_new::new;
use itertools::Itertools;
use ordered_float::NotNan;
use rust_htslib::bcf::{self, record::GenotypeAllele::Unphased, Read};
use std::collections::BTreeMap;
use std::iter::FromIterator;
use std::{path::PathBuf, str};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    haplotype_variants: bcf::Reader,
    haplotype_calls: bcf::Reader,
    observations: bcf::Reader,
    max_haplotypes: usize,
    outcsv: Option<PathBuf>,
    k_reads: usize,
}

impl Caller {
    pub fn call(&mut self) -> Result<()> {
        // Step 1: obtain estimates and varlociraptor.
        //1) loop through the different HLA loci
        //for testing, DQA1 haplotypes
        // let haplotypes = vec![
        //     "HLA:HLA21879",
        //     "HLA:HLA06599",
        //     "HLA:HLA00604",
        //     "HLA:HLA00612",
        //     "HLA:HLA00617",
        //     "HLA:HLA00610",
        //     "HLA:HLA02433",
        //     "HLA:HLA01409",
        //     "HLA:HLA14846",
        //     "HLA:HLA06601",
        //     "HLA:HLA06598",
        //     "HLA:HLA14797",
        //     "HLA:HLA01376",
        //     "HLA:HLA00602",
        //     "HLA:HLA00611",
        //     "HLA:HLA21186",
        //     "HLA:HLA06618",
        //     "HLA:HLA00606",
        //     "HLA:HLA00608",
        //     "HLA:HLA00603",
        //     "HLA:HLA01694",
        //     "HLA:HLA02339",
        //     "HLA:HLA00607",
        //     "HLA:HLA17309",
        //     "HLA:HLA00601",
        //     "HLA:HLA01905",
        //     "HLA:HLA00605",
        //     "HLA:HLA00613",
        //     "HLA:HLA00619",
        //     "HLA:HLA06614",
        //     "HLA:HLA02340",
        //     "HLA:HLA09875",
        //     "HLA:HLA19774",
        //     "HLA:HLA01697",
        //     "HLA:HLA06609",
        // ];
        // let haplotypes: Vec<Haplotype> = haplotypes
        //     .iter()
        //     .map(|haplotype| Haplotype(haplotype.to_string()))
        //     .collect();
        //1) obtain varlociraptor calls according to criteria:
        //read_depths > 0, afd field is not empty and prob_absent <= 0.1.
        let mut haplotype_calls = HaplotypeCalls::new(&mut self.haplotype_calls)?;

        //collect the variant IDs from the varlociraptor calls.
        let variant_ids: Vec<VariantID> = haplotype_calls.keys().cloned().collect();
        //2) collect the candidate variants and haplotypes (horizontal evidence, shortlist of haplotypes)
        let haplotype_variants = HaplotypeVariants::new(
            &mut self.observations,
            &mut self.haplotype_variants,
            &variant_ids,
            &self.k_reads,
            &self.max_haplotypes,
        )?;
        //4)keep only the variants secondly filtered in haplotype_variants
        haplotype_calls.retain(|&k, _| haplotype_variants.contains_key(&k));

        //5) finalize haplotype variants struct to contain only thee shortlisted haplotypes and
        //select top N haplotypes according to --max-haplotypes N.
        let (_, haplotypes_gt_c) = haplotype_variants.iter().next().unwrap();
        let haplotypes: Vec<Haplotype> = haplotypes_gt_c.keys().cloned().collect();
        //6) create the GenotypesLoci struct that contains variants, genotypes and loci information,
        //together with only selected top N haplotypes from the previous step.
        let variant_matrix = GenotypesLoci::new(&haplotype_variants, &haplotypes).unwrap();
        // Step 2: setup model.
        let model = Model::new(Likelihood::new(), Prior::new(), Posterior::new());
        let data = Data::new(variant_matrix, haplotype_calls);
        // Step 3: calculate posteriors.
        //let m = model.compute(universe, &data);
        let m = model.compute_from_marginal(&Marginal::new(haplotypes.len()), &data);
        let posterior_output = m.event_posteriors();

        //add variant query and probabilities to the outout table for each event
        let variant_calls: Vec<AlleleFreqDist> = data.haplotype_calls.values().cloned().collect();
        let genotype_loci_matrix = data.variant_matrix;
        let mut event_queries: Vec<BTreeMap<VariantID, (AlleleFreq, LogProb)>> = Vec::new();
        posterior_output.for_each(|(fractions, _)| {
            let mut vaf_queries: BTreeMap<VariantID, (AlleleFreq, LogProb)> = BTreeMap::new();
            genotype_loci_matrix
                .iter()
                .zip(variant_calls.iter())
                .for_each(|((variant_id, (genotypes, covered)), afd)| {
                    let mut denom = NotNan::new(1.0).unwrap();
                    let mut vaf_sum = NotNan::new(0.0).unwrap();
                    let mut counter = 0;
                    fractions.iter().enumerate().for_each(|(i, fraction)| {
                        if genotypes[i as u64] && covered[i as u64] {
                            vaf_sum += *fraction;
                            counter += 1;
                        } else if covered[i as u64] {
                            ()
                        } else {
                            denom -= *fraction;
                        }
                    });
                    if denom > NotNan::new(0.0).unwrap() {
                        vaf_sum /= denom;
                    }
                    vaf_sum = NotNan::new((vaf_sum * NotNan::new(100.0).unwrap()).round()).unwrap()
                        / NotNan::new(100.0).unwrap();
                    let answer = afd.vaf_query(&vaf_sum);
                    if counter > 0 {
                        vaf_queries.insert(*variant_id, (vaf_sum, answer));
                    }
                });
            event_queries.push(vaf_queries);
        });
        // Step 4: print TSV table with results
        // TODO use csv crate
        // Columns: posterior_prob, haplotype_a, haplotype_b, haplotype_c, ...
        // with each column after the first showing the fraction of the respective haplotype
        let mut wtr = csv::Writer::from_path(self.outcsv.as_ref().unwrap())?;
        let mut headers: Vec<_> = vec!["density".to_string(), "odds".to_string()];
        let haplotypes_str: Vec<String> = haplotypes.iter().map(|h| h.to_string()).collect();
        headers.extend(haplotypes_str);
        let variant_names = event_queries[0]
            .keys()
            .map(|key| format!("{:?}", key))
            .collect::<Vec<String>>();
        headers.extend(variant_names); //add variant names as separate columns
        wtr.write_record(&headers)?;

        //write best record on top
        let mut records = Vec::new();
        let mut posterior = m.event_posteriors();
        let (haplotype_frequencies, best_density) = posterior.next().unwrap();
        let best_odds = 1;
        let format_f64 = |number: f64, records: &mut Vec<String>| {
            if number <= 0.01 {
                records.push(format!("{:+.2e}", number))
            } else {
                records.push(format!("{:.2}", number))
            }
        };
        format_f64(best_density.exp(), &mut records);
        records.push(best_odds.to_string());
        let format_freqs = |frequency: NotNan<f64>, records: &mut Vec<String>| {
            if frequency <= NotNan::new(0.01).unwrap() {
                records.push(format!("{:+.2e}", NotNan::into_inner(frequency)))
            } else {
                records.push(format!("{:.2}", frequency))
            }
        };
        haplotype_frequencies
            .iter()
            .for_each(|frequency| format_freqs(*frequency, &mut records));
        //add vaf queries and probabilities for the first event to the output table
        let queries: Vec<(AlleleFreq, LogProb)> = event_queries
            .iter()
            .next()
            .unwrap()
            .values()
            .cloned()
            .collect();
        queries.iter().for_each(|(query, answer)| {
            let prob = f64::from(Prob::from(*answer));
            if prob <= 0.01 {
                records.push(format!("{}{}{:+.2e}", query, ":", prob));
            } else {
                records.push(format!("{}{}{:.2}", query, ":", prob));
            }
        });
        wtr.write_record(records)?;

        //write the rest of the records
        posterior.zip(event_queries.iter().skip(1)).for_each(
            |((haplotype_frequencies, density), queries)| {
                let mut records = Vec::new();
                let odds = (density - best_density).exp();
                format_f64(density.exp(), &mut records);
                format_f64(odds, &mut records);
                haplotype_frequencies
                    .iter()
                    .for_each(|frequency| format_freqs(*frequency, &mut records));

                queries.iter().for_each(|(_, (query, answer))| {
                    let prob = f64::from(Prob::from(*answer));
                    if prob <= 0.01 {
                        records.push(format!("{}{}{:+.2e}", query, ":", prob));
                    } else {
                        records.push(format!("{}{}{:.2}", query, ":", prob));
                    }
                });
                wtr.write_record(records).unwrap();
            },
        );
        Ok(())
    }
}

#[derive(Derefable, Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, new)]
pub(crate) struct Haplotype(#[deref] String);

#[derive(Derefable, Debug, Copy, Clone, PartialEq, Eq, Hash, Ord, PartialOrd, new)]
pub(crate) struct VariantID(#[deref] i32);

#[derive(Derefable, Debug, Clone, PartialEq, Eq, PartialOrd)]
pub(crate) struct HaplotypeVariants(
    #[deref] BTreeMap<VariantID, BTreeMap<Haplotype, (bool, bool)>>,
);

impl HaplotypeVariants {
    pub(crate) fn new(
        observations: &mut bcf::Reader,
        haplotype_variants: &mut bcf::Reader,
        filtered_ids: &Vec<VariantID>,
        k_reads: &usize,
        max_haplotypes: &usize,
    ) -> Result<Self> {
        let mut variant_fragment_map: BTreeMap<bool, Vec<u64>> = BTreeMap::new();
        let mut variant_records = BTreeMap::new(); //these two maps have the same variants in the same order so only variant_records keeps their name to prevent redundant disk usage.
                                                   //collection of fragments
        for record_result in observations.records() {
            let mut record = record_result?;
            let variant_id: VariantID = VariantID(String::from_utf8(record.id())?.parse().unwrap());
            if filtered_ids.contains(&variant_id) {
                let read_observations =
                    varlociraptor::calling::variants::preprocessing::read_observations(&mut record)
                        .unwrap();
                let read_observation = &read_observations.pileup.read_observations();
                read_observation.iter().for_each(|read| {
                    let fragment_id = read.fragment_id.unwrap();
                    let alt_evidence = read.prob_alt > read.prob_ref;
                    variant_fragment_map.entry(alt_evidence).or_insert(vec![]);
                    let mut fragments = variant_fragment_map.get(&alt_evidence).unwrap().clone();
                    fragments.push(fragment_id);
                    variant_fragment_map.insert(alt_evidence, fragments);
                })
            }
        }
        dbg!(&variant_fragment_map);
        //collection of variant matrices
        for record_result in haplotype_variants.records() {
            let record = record_result?;
            let variant_id: VariantID = VariantID(String::from_utf8(record.id())?.parse().unwrap());
            if filtered_ids.contains(&variant_id) {
                let header = record.header();
                let gts = record.genotypes()?;
                let loci = record.format(b"C").integer().unwrap();
                let mut matrices = BTreeMap::new();
                for (index, haplotype) in header.samples().iter().enumerate() {
                    let haplotype = Haplotype(str::from_utf8(haplotype).unwrap().to_string());
                    for gta in gts.get(index).iter() {
                        matrices.insert(
                            haplotype.clone(),
                            ((*gta == Unphased(1)), (loci[index] == &[1])),
                        );
                    }
                }
                variant_records.insert(variant_id, matrices);
            }
        }
        dbg!(&variant_records);
        //loop over variant fragment map and variant records at the same time.
        let mut allele_count: BTreeMap<Haplotype, u64> = BTreeMap::new();
        variant_fragment_map
            .iter()
            .zip(variant_records.iter())
            .for_each(|((alt_evidence, fragments), (variant, matrices))| {
                dbg!(&variant);
                fragments.iter().for_each(|fragment| {
                    let mut variant_in_haplotype = false;
                    let mut last_visited_haplotype = String::from("");
                    for (haplotype, (genotype, locus)) in matrices {
                        if *genotype {
                            let chk2 = 1;
                            dbg!(&chk2);
                            dbg!(&variant);
                            allele_count.entry(haplotype.clone()).or_insert(0);
                            let mut count = allele_count.get(&haplotype).unwrap().clone();
                            count += 1;
                            if variant_in_haplotype {
                                let mut previous_haplotype_count = allele_count
                                    .get(&Haplotype(last_visited_haplotype.clone()))
                                    .unwrap()
                                    .clone();
                                previous_haplotype_count -= 1;
                                allele_count.insert(
                                    Haplotype(last_visited_haplotype),
                                    previous_haplotype_count,
                                );
                                break;
                            } else {
                                allele_count.insert(haplotype.clone(), count);
                                last_visited_haplotype.push_str(&haplotype.to_string()); //push the name to the memory as the last visited haplotype
                                variant_in_haplotype = true;
                            }
                        }
                    }
                });
            });
        dbg!(&allele_count);
        //sort the map and get --max-n-haplotypes
        let mut v = Vec::from_iter(allele_count);
        v.sort_by(|(_, a), (_, b)| b.cmp(&a));
        dbg!(&v);
        let max_n_haplotypes: Vec<Haplotype> = v[0..*max_haplotypes]
            .iter()
            .map(|(haplotype, _)| haplotype.clone())
            .collect();

        // //create two maps: 1) a variantID, fragments map and a variantID, haplotypes map.
        // let mut variants_fragments: BTreeMap<VariantID, Vec<u64>> = BTreeMap::new();
        // let mut variants_haplotypes: BTreeMap<VariantID, Vec<Haplotype>> = BTreeMap::new();
        // let mut variant_records = BTreeMap::new();
        // for record_result in observations.records() {
        //     let mut record = record_result?;
        //     let variant_id: VariantID = VariantID(String::from_utf8(record.id())?.parse().unwrap());
        //     if filtered_ids.contains(&variant_id) {
        //         let read_observations =
        //             varlociraptor::calling::variants::preprocessing::read_observations(&mut record)
        //                 .unwrap();
        //         let read_observation = &read_observations.pileup.read_observations();
        //         read_observation.iter().for_each(|read| {
        //             if read.prob_ref > read.prob_alt {
        //                 //read.prob_alt > read.prob_ref
        //                 variants_fragments.entry(variant_id).or_insert(vec![]);
        //                 let mut fragments = variants_fragments.get(&variant_id).unwrap().clone();
        //                 fragments.push(read.fragment_id.unwrap());
        //                 variants_fragments.insert(variant_id, fragments); //overwrite the existing key-value pair with the key-value update.
        //             }
        //         })
        //     }
        // }
        // for record_result in haplotype_variants.records() {
        //     let record = record_result?;
        //     let variant_id: VariantID = VariantID(String::from_utf8(record.id())?.parse().unwrap());
        //     if filtered_ids.contains(&variant_id) {
        //         let header = record.header();
        //         let gts = record.genotypes()?;
        //         let loci = record.format(b"C").integer().unwrap();
        //         let mut matrices = BTreeMap::new();

        //         for (index, haplotype) in header.samples().iter().enumerate() {
        //             let haplotype = Haplotype(str::from_utf8(haplotype).unwrap().to_string());
        //             if haplotypes.contains(&haplotype) {
        //                 for gta in gts.get(index).iter() {
        //                     matrices.insert(
        //                         haplotype.clone(),
        //                         ((*gta == Unphased(1)), (loci[index] == &[1])),
        //                     );
        //                     if loci[index] == &[1] {
        //                         //*gta == Unphased(1)
        //                         variants_haplotypes.entry(variant_id).or_insert(vec![]);
        //                         let mut haplotypes =
        //                             variants_haplotypes.get(&variant_id).unwrap().clone();
        //                         haplotypes.push(haplotype.clone());
        //                         variants_haplotypes.insert(variant_id, haplotypes);
        //                     }
        //                 }
        //             }
        //         }
        //         variant_records.insert(variant_id, matrices);
        //     }
        // }
        // //count the number of variants each fragment bears.
        // let mut fragment_count: BTreeMap<u64, u64> = BTreeMap::new();
        // variants_fragments.iter().for_each(|(_, fragments)| {
        //     fragments.iter().for_each(|fragment| {
        //         fragment_count.entry(*fragment).or_insert(0);
        //         let mut count = fragment_count.get(&fragment).unwrap().clone();
        //         count += 1;
        //         fragment_count.insert(*fragment, count);
        //     });
        // });
        // let mut filtered_haplotypes: Vec<Haplotype> = Vec::new();
        // variants_haplotypes
        //     .iter()
        //     .for_each(|(variant, haplotypes)| {
        //         if let Some(fragments) = variants_fragments.get(&variant) {
        //             if fragments.len() > *k_reads {
        //                 if fragments
        //                     .iter()
        //                     .all(|fragment| fragment_count.get(&fragment).unwrap() >= &2)
        //                 {
        //                     filtered_haplotypes.extend(haplotypes.clone());
        //                 }
        //             }
        //         }
        //     });
        // let final_haplotypes: Vec<Haplotype> = filtered_haplotypes.into_iter().unique().collect();
        // //count the number of variants each haplotype bears.
        // let mut haplotype_count: BTreeMap<Haplotype, u64> = BTreeMap::new();
        // variants_haplotypes
        //     .iter()
        //     .for_each(|(variant, haplotypes)| {
        //         haplotypes.iter().for_each(|haplotype| {
        //             haplotype_count.entry(haplotype.clone()).or_insert(0);
        //             let mut count = haplotype_count.get(&haplotype).unwrap().clone();
        //             count += 1;
        //             haplotype_count.insert(haplotype.clone(), count);
        //         });
        //     });
        // //subset only the haplotypes that are contained in final_haplotypes
        // let filtered_haplotype_count: BTreeMap<Haplotype, u64> = haplotype_count
        //     .iter()
        //     .filter(|(haplotype, _)| final_haplotypes.contains(&haplotype.clone()))
        //     .map(|(haplotype, count)| (haplotype.clone(), *count))
        //     .collect();
        //sort the map and get --max-n-haplotypes
        // let mut v = Vec::from_iter(filtered_haplotype_count);
        // v.sort_by(|(_, a), (_, b)| b.cmp(&a));
        // //dbg!(&v);
        // let max_n_haplotypes: Vec<Haplotype> = v[0..*max_haplotypes]
        //     .iter()
        //     .map(|(haplotype, _)| haplotype.clone())
        //     .collect();

        //keep only the variants of variant_fragments (bc of prob_alt>prob_ref condition)
        //variant_records.retain(|&k, _| variants_fragments.contains_key(&k));
        //filter both for the selected haplotypes and the variants.
        let mut variant_records = HaplotypeVariants(variant_records);
        let filtered_variant_records =
            HaplotypeVariants::filter_haplotypes(&mut variant_records, &max_n_haplotypes).unwrap();
        Ok(filtered_variant_records)
    }
    fn filter_haplotypes(
        haplotype_variants: &mut Self,
        filtered_haps: &Vec<Haplotype>,
    ) -> Result<Self> {
        let mut filtered_haplotypes: Vec<String> = Vec::new();
        //the loop for discovering haplotype names that bear the filtered variants.
        for (_, genotypes_loci_map) in haplotype_variants.iter() {
            for (haplotype, (genotype, _)) in genotypes_loci_map {
                if *genotype && filtered_haps.contains(haplotype) {
                    filtered_haplotypes.push(haplotype.to_string());
                }
            }
        }
        //remove the duplicates.
        let filtered_haplotypes: Vec<String> = filtered_haplotypes.into_iter().unique().collect();
        //1) collect the first record of haplotype_variants and collect the indices of haplotypes for further filtering in the following.
        let mut haplotype_indices: Vec<usize> = Vec::new();
        if let Some((_, bmap)) = haplotype_variants.iter().next() {
            bmap.iter().enumerate().for_each(|(i, (haplotype, _))| {
                if filtered_haplotypes.contains(haplotype) {
                    haplotype_indices.push(i)
                }
            });
        };
        //create the filtered matrix with the indices of filtered haplotypes.
        let mut new_variants = BTreeMap::new();
        haplotype_variants.iter().for_each(|(variant_id, bmap)| {
            let mut matrix_map = BTreeMap::new();
            bmap.iter()
                .enumerate()
                .for_each(|(i, (haplotype, matrices))| {
                    haplotype_indices.iter().for_each(|j| {
                        if i == *j {
                            matrix_map.insert(haplotype.clone(), *matrices);
                        }
                    });
                });
            new_variants.insert(*variant_id, matrix_map);
        });
        Ok(HaplotypeVariants(new_variants))
    }
}
#[derive(Debug, Clone, Derefable)]
pub(crate) struct AlleleFreqDist(#[deref] BTreeMap<AlleleFreq, f64>);

impl AlleleFreqDist {
    pub(crate) fn vaf_query(&self, vaf: &AlleleFreq) -> LogProb {
        if self.contains_key(&vaf) {
            LogProb::from(PHREDProb(*self.get(&vaf).unwrap()))
        } else {
            let (x_0, y_0) = self.range(..vaf).next_back().unwrap();
            let (x_1, y_1) = self.range(vaf..).next().unwrap();
            let density =
                NotNan::new(*y_0).unwrap() + (*vaf - *x_0) * (*y_1 - *y_0) / (*x_1 - *x_0); //calculation of density for given vaf by linear interpolation
            LogProb::from(PHREDProb(NotNan::into_inner(density)))
        }
    }
}

#[derive(Derefable, Debug)]
pub(crate) struct GenotypesLoci(#[deref] BTreeMap<VariantID, (BitVec, BitVec)>);

impl GenotypesLoci {
    pub(crate) fn new(
        haplotype_variants: &HaplotypeVariants,
        haplotypes: &Vec<Haplotype>,
    ) -> Result<Self> {
        let mut variant_matrix = BTreeMap::new();
        haplotype_variants.iter().for_each(|(variant_id, bmap)| {
            let mut haplotype_variants_gt = BitVec::new();
            let mut haplotype_variants_c = BitVec::new();
            bmap.iter().for_each(|(haplotype, (gt, c))| {
                haplotypes.iter().for_each(|h| {
                    if haplotype == h {
                        haplotype_variants_gt.push(*gt);
                        haplotype_variants_c.push(*c);
                    }
                });
            });
            variant_matrix.insert(*variant_id, (haplotype_variants_gt, haplotype_variants_c));
        });
        Ok(GenotypesLoci(variant_matrix))
    }
}

#[derive(Derefable, DerefMut, Debug, Clone)]
pub(crate) struct HaplotypeCalls(#[deref] BTreeMap<VariantID, AlleleFreqDist>);

impl HaplotypeCalls {
    pub(crate) fn new(haplotype_calls: &mut bcf::Reader) -> Result<Self> {
        let mut calls = BTreeMap::new();
        for record_result in haplotype_calls.records() {
            let record = record_result?;
            let prob_absent = record.info(b"PROB_ABSENT").float().unwrap().unwrap()[0];
            let prob_absent_prob = Prob::from(PHREDProb(prob_absent.into()));
            let afd_utf = record.format(b"AFD").string()?;
            let afd = std::str::from_utf8(afd_utf[0]).unwrap();
            let read_depths = record.format(b"DP").integer().unwrap();
            let variant_id: i32 = String::from_utf8(record.id())?.parse().unwrap();
            if read_depths[0] != &[0] && afd != "." && &prob_absent_prob <= &Prob(0.1)
                || &prob_absent_prob >= &Prob(0.9)
            {
                //because some afd strings are just "." and that throws an error while splitting below.
                let variant_id: i32 = String::from_utf8(record.id())?.parse().unwrap();
                let mut vaf_density = BTreeMap::new();
                for pair in afd.split(',') {
                    let (vaf, density) = pair.split_once("=").unwrap();
                    let (vaf, density): (AlleleFreq, f64) =
                        (vaf.parse().unwrap(), density.parse().unwrap());
                    vaf_density.insert(vaf, density);
                }
                calls.insert(VariantID(variant_id), AlleleFreqDist(vaf_density));
            }
        }
        Ok(HaplotypeCalls(calls))
    }
}
