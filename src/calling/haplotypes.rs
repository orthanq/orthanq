use crate::model::{AlleleFreq, Data, Likelihood, Marginal, Posterior, Prior};
use anyhow::Result;
use bio::stats::{bayesian::model::Model, probs::LogProb, PHREDProb, Prob};
use bv::BitVec;
use derefable::Derefable;
use derive_builder::Builder;
use derive_deref::DerefMut;
use derive_new::new;
use hdf5;
use itertools::Itertools;
use ordered_float::NotNan;
use rust_htslib::bcf::{self, record::GenotypeAllele::Unphased, Read};
use std::collections::BTreeMap;
use std::convert::TryInto;
use std::{path::PathBuf, str};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    hdf5_reader: hdf5::File,
    haplotype_variants: bcf::Reader,
    haplotype_calls: bcf::Reader,
    min_norm_counts: f64,
    max_haplotypes: i64,
    outcsv: Option<PathBuf>,
    use_evidence: String,
}

impl Caller {
    pub fn call(&mut self) -> Result<()> {
        // Step 1: obtain kallisto estimates.
        //obtain the varlociraptor calls.
        let haplotype_calls = HaplotypeCalls::new(&mut self.haplotype_calls)?;
        dbg!(&haplotype_calls.len());
        //extract the filtered variant IDs according to read_depths > 0, afd field is not empty and prob_absent <= 0.1.
        let filtered_ids: Vec<VariantID> = haplotype_calls.keys().cloned().collect();
        dbg!(&filtered_ids.len());
        //filter the candidate variants and retain only the selected haplotypes.
        let haplotype_variants =
            HaplotypeVariants::new(&mut self.haplotype_variants, &filtered_ids)?;
        dbg!(&haplotype_variants.len());
        //collect the haplotype names for the creation of KallistoEstimates with only the selected haplotypes.
        let vals: Vec<BTreeMap<String, (bool, bool)>> =
            haplotype_variants.values().cloned().collect();
        let haplotypes: Vec<String> = vals[0].keys().cloned().collect();

        //create the KallistoEstimates struct with pre-filtered haplotypes.
        let kallisto_estimates = KallistoEstimates::new(
            &self.hdf5_reader,
            self.min_norm_counts,
            self.max_haplotypes,
            &haplotypes,
        )?;
        dbg!(&kallisto_estimates);
        //collect the names of final N number of haplotypes according to --max-haplotypes
        let final_haplotypes: Vec<String> =
            kallisto_estimates.keys().map(|x| x.to_string()).collect();

        //Create the GenotypesLoci struct that contains variants, genotypes and loci information.
        let variant_matrix = GenotypesLoci::new(&haplotype_variants, &final_haplotypes).unwrap();
        dbg!(&variant_matrix.len());
        // Step 2: setup model.
        let model = Model::new(
            Likelihood::new(self.use_evidence.clone()),
            Prior::new(),
            Posterior::new(),
        );
        let data = Data::new(
            kallisto_estimates.values().cloned().collect(),
            variant_matrix,
            haplotype_calls,
        );
        // Step 3: calculate posteriors.
        //let m = model.compute(universe, &data);
        dbg!(&final_haplotypes.len());
        let m = model.compute_from_marginal(&Marginal::new(final_haplotypes.len()), &data);
        let posterior_output = m.event_posteriors();
        //add variant query and probabilities to the outout table for each event
        //let variant_matrix: Vec<(BitVec, BitVec)> = data.variant_matrix.values().cloned().collect();
        //dbg!(&variant_matrix);
        let variant_calls: Vec<AlleleFreqDist> = data.haplotype_calls.values().cloned().collect();
        let genotype_loci_matrix = data.variant_matrix;
        let mut event_queries: Vec<BTreeMap<VariantID, (AlleleFreq, LogProb)>> = Vec::new();
        posterior_output.for_each(|(fractions, _)| {
            let mut vaf_queries: BTreeMap<VariantID, (AlleleFreq, LogProb)> = BTreeMap::new();
            genotype_loci_matrix.iter().zip(variant_calls.iter()).for_each(
                |((variant_id,(genotypes, covered)), afd)| {
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
                    dbg!(&variant_id);
                    dbg!(&vaf_sum);
                    dbg!(&answer);
                },
            );
            event_queries.push(vaf_queries);
        });
        // Step 4: print TSV table with results
        // TODO use csv crate
        // Columns: posterior_prob, haplotype_a, haplotype_b, haplotype_c, ...
        // with each column after the first showing the fraction of the respective haplotype
        let mut wtr = csv::Writer::from_path(self.outcsv.as_ref().unwrap())?;
        let mut headers: Vec<_> = vec!["density".to_string(), "odds".to_string()];
        headers.extend(final_haplotypes);
        let variant_names = event_queries[0]
            .keys()
            .map(|key| format!("{:?}",key))
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

#[derive(Derefable, Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) struct Haplotype(#[deref] String);

#[derive(Debug, Clone, Copy, Ord, PartialOrd, Eq, PartialEq)]
pub(crate) struct KallistoEstimate {
    pub count: NotNan<f64>,
    pub dispersion: NotNan<f64>,
}

#[derive(Debug, Clone, Derefable)]
pub(crate) struct KallistoEstimates(#[deref] BTreeMap<Haplotype, KallistoEstimate>);

impl KallistoEstimates {
    /// Generate new instance.
    pub(crate) fn new(
        hdf5_reader: &hdf5::File,
        min_norm_counts: f64,
        max_haplotypes: i64,
        haplotypes: &Vec<String>,
    ) -> Result<Self> {
        let seqnames = Self::filter_seqnames(hdf5_reader, min_norm_counts, &haplotypes)?;
        let ids = hdf5_reader
            .dataset("aux/ids")?
            .read_1d::<hdf5::types::FixedAscii<255>>()?;
        let num_bootstraps = hdf5_reader.dataset("aux/num_bootstrap")?.read_1d::<i32>()?;
        let seq_length = hdf5_reader.dataset("aux/lengths")?.read_1d::<f64>()?;
        let mut estimates = BTreeMap::new();
        for seqname in seqnames {
            let index = ids.iter().position(|x| x.as_str() == seqname).unwrap();
            let mut bootstraps = Vec::new();
            for i in 0..num_bootstraps[0] {
                let dataset = hdf5_reader.dataset(&format!("bootstrap/bs{i}", i = i))?;
                let est_counts = dataset.read_1d::<f64>()?;
                let norm_counts = est_counts / &seq_length;
                let norm_counts = norm_counts[index];
                bootstraps.push(norm_counts);
            }

            //mean
            let sum = bootstraps.iter().sum::<f64>();
            let count = bootstraps.len();
            let m = sum / count as f64;

            //std dev
            let variance = bootstraps
                .iter()
                .map(|value| {
                    let diff = m - (*value as f64);
                    diff * diff
                })
                .sum::<f64>()
                / count as f64;
            let std = variance.sqrt();
            let t = std / m;
            dbg!(&seqname);
            dbg!(&std);
            dbg!(&m);
            dbg!(&t);
            //retrieval of mle
            let mle_dataset = hdf5_reader.dataset("est_counts")?.read_1d::<f64>()?;
            let mle_norm = mle_dataset / &seq_length; //normalized mle counts by length
            let m = mle_norm[index];
            estimates.insert(
                Haplotype(seqname.clone()),
                KallistoEstimate {
                    dispersion: NotNan::new(t).unwrap(),
                    count: NotNan::new(m).unwrap(),
                },
            );
        }
        let kallisto_estimates = KallistoEstimates(estimates);
        Self::select_haplotypes(kallisto_estimates, max_haplotypes)
    }

    //Return top N estimates according to --max-haplotypes
    fn select_haplotypes(self, max_haplotypes: i64) -> Result<Self> {
        let mut estimates_vec: Vec<(&Haplotype, &KallistoEstimate)> = self.iter().collect();
        estimates_vec.sort_by(|a, b| b.1.count.cmp(&a.1.count));
        if estimates_vec.len() >= max_haplotypes.try_into().unwrap() {
            let topn = estimates_vec[0..max_haplotypes as usize].to_vec();
            let mut top_estimates = BTreeMap::new();
            for (key, value) in topn {
                top_estimates.insert(key.clone(), *value);
            }
            Ok(KallistoEstimates(top_estimates))
        } else {
            Ok(self)
        }
    }
    //Return a vector of filtered seqnames according to --min-norm-counts.
    fn filter_seqnames(
        hdf5_reader: &hdf5::File,
        min_norm_counts: f64,
        haplotypes: &Vec<String>,
    ) -> Result<Vec<String>> {
        let ids = hdf5_reader
            .dataset("aux/ids")?
            .read_1d::<hdf5::types::FixedAscii<255>>()?;
        let mut haplotype_indices = Vec::new();
        ids.iter().enumerate().for_each(|(index, id)| {
            if haplotypes.contains(&id.as_str().to_string()) {
                haplotype_indices.push(index as usize);
            }
        });
        let est_counts = hdf5_reader.dataset("est_counts")?.read_1d::<f64>()?;
        let seq_length = hdf5_reader.dataset("aux/lengths")?.read_1d::<f64>()?; //these two variables arrays have the same length.
        let norm_counts = est_counts / seq_length;
        let mut norm_counts_filtered = Vec::new();
        haplotype_indices
            .iter()
            .for_each(|i| norm_counts_filtered.push(norm_counts[*i]));
        let mut filtered_indices = Vec::new();
        for (num, index) in norm_counts_filtered.iter().zip(haplotype_indices) {
            if num > &min_norm_counts {
                filtered_indices.push(index);
            }
        }
        let mut filtered_haplotypes: Vec<String> = Vec::new();
        for i in filtered_indices {
            filtered_haplotypes.push(ids[i].to_string());
        }
        Ok(filtered_haplotypes)
    }
}

#[derive(Derefable, Debug, Copy, Clone, PartialEq, Eq, Hash, Ord, PartialOrd, new)]
pub(crate) struct VariantID(#[deref] i32);

#[derive(Derefable, DerefMut, Debug, Clone, PartialEq, Eq, PartialOrd)]
pub(crate) struct HaplotypeVariants(#[deref] BTreeMap<VariantID, BTreeMap<String, (bool, bool)>>);

impl HaplotypeVariants {
    pub(crate) fn new(
        haplotype_variants: &mut bcf::Reader,
        filtered_ids: &Vec<VariantID>,
    ) -> Result<Self> {
        let mut variant_records = BTreeMap::new();
        for record_result in haplotype_variants.records() {
            let record = record_result?;
            let variant_id: i32 = String::from_utf8(record.id())?.parse().unwrap();
            if filtered_ids.contains(&VariantID(variant_id)) {
                let header = record.header();
                let gts = record.genotypes()?; //genotypes of all samples
                let loci = record.format(b"C").integer().unwrap();
                let mut matrices = BTreeMap::new();
                //let mut all_haplotypes: Vec<String> = Vec::new();
                for (index, haplotype) in header.samples().iter().enumerate() {
                    let haplotype_name = str::from_utf8(haplotype).unwrap().to_string();
                    for gta in gts.get(index).iter() {
                        //dbg!(&gta);
                        //dbg!(&locus);
                        matrices.insert(
                            haplotype_name.clone(),
                            ((*gta == Unphased(1)), (loci[index] == &[1])),
                        );
                    }
                }
                variant_records.insert(VariantID(variant_id), matrices);
            }
        }
        let mut variant_records = HaplotypeVariants(variant_records);
        //filter both for the selected haplotypes and the variants.
        HaplotypeVariants::filter_haplotypes(&mut variant_records, &filtered_ids).unwrap();
        Ok(variant_records)
    }
    fn filter_haplotypes(
        haplotype_variants: &mut Self,
        filtered_variants: &Vec<VariantID>,
    ) -> Result<Self> {
        let mut filtered_haplotypes: Vec<String> = Vec::new();
        //the loop for discovering haplotype names that bear the filtered variants.
        for (variant, genotypes_loci_map) in haplotype_variants.iter() {
            for (haplotype, (genotype, _)) in genotypes_loci_map {
                if filtered_variants.contains(variant) && *genotype {
                    filtered_haplotypes.push(haplotype.to_string());
                }
            }
        }
        //remove the duplicates.
        let filtered_haplotypes: Vec<String> = filtered_haplotypes.into_iter().unique().collect();
        dbg!(&filtered_haplotypes.len());
        //collect indices of filtered_haplotypes.
        let mut haplotype_indices: Vec<usize> = Vec::new();
        //1) collect the first record of haplotype_variants.
        let mut bmaptree = BTreeMap::new();
        for (_, bmap) in haplotype_variants.iter() {
            for (k, v) in bmap.iter() {
                bmaptree.insert(k, v);
            }
            break;
        }
        dbg!(&haplotype_variants.len());
        //2) then collect the indices of haplotypes for further filtering in the following.
        bmaptree.iter().enumerate().for_each(|(i, (haplotype, _))| {
            if filtered_haplotypes.contains(haplotype) {
                haplotype_indices.push(i)
            }
        });
        dbg!(&haplotype_indices.len());
        //create the filtered matrix with the indices of filtered haplotypes.
        let mut new_variants = BTreeMap::new();
        haplotype_variants.iter().for_each(|(variant_id, bmap)| {
            if filtered_variants.contains(variant_id) {
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
            }
        });
        dbg!(&new_variants.len());
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
        haplotypes: &Vec<String>,
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
        dbg!(&variant_matrix);
        Ok(GenotypesLoci(variant_matrix))
    }
}

#[derive(Derefable, Debug, Clone)]
pub(crate) struct HaplotypeCalls(#[deref] BTreeMap<VariantID, AlleleFreqDist>);

impl HaplotypeCalls {
    pub(crate) fn new(haplotype_calls: &mut bcf::Reader) -> Result<Self> {
        let mut calls = BTreeMap::new();
        for record_result in haplotype_calls.records() {
            let record = record_result?;
            let prob_absent = record.info(b"PROB_ABSENT").float().unwrap().unwrap()[0];
            let variant_id: i32 = String::from_utf8(record.id())?.parse().unwrap();
            dbg!(&variant_id);
            dbg!(&prob_absent);
            let prob_absent_prob = Prob::from(PHREDProb(prob_absent.into()));
            dbg!(&prob_absent_prob);
            let afd_utf = record.format(b"AFD").string()?;
            let afd = std::str::from_utf8(afd_utf[0]).unwrap();
            let read_depths = record.format(b"DP").integer().unwrap();
            if read_depths[0] != &[0] && afd != "." && &prob_absent_prob <= &Prob(0.1) {
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
