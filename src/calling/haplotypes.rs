use crate::model::{AlleleFreq, Data, HaplotypeFractions, Likelihood, Marginal, Posterior, Prior};
use anyhow::Result;
use bio::stats::{bayesian::model::Model, probs::LogProb, PHREDProb, Prob};
use bv::BitVec;
use derefable::Derefable;
use derive_builder::Builder;
use derive_deref::DerefMut;
use derive_new::new;
use itertools::Itertools;
use linfa::prelude::*;
use linfa_clustering::KMeans;
use ndarray::prelude::*;
use ordered_float::NotNan;
use plotters::prelude::*;
use rand::prelude::*;
use rust_htslib::bcf::{self, record::GenotypeAllele::Unphased, Read};
use std::collections::BTreeMap;
use std::{path::PathBuf, str};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    haplotype_variants: bcf::Reader,
    variant_calls: bcf::Reader,
    max_haplotypes: usize,
    outcsv: Option<PathBuf>,
}

impl Caller {
    pub fn call(&mut self) -> Result<()> {
        let mut variant_calls = VariantCalls::new(&mut self.variant_calls)?;
        let variant_ids: Vec<VariantID> = variant_calls.keys().cloned().collect();
        let mut haplotype_variants =
            HaplotypeVariants::new(&mut self.haplotype_variants, &variant_ids)?;
        let haplotype_variants =
            haplotype_variants.find_plausible_haplotypes(self.max_haplotypes)?;
        let (_, candidate_matrix) = haplotype_variants.iter().next().unwrap();
        let haplotypes: Vec<Haplotype> = candidate_matrix.keys().cloned().collect();
        let candidate_matrix = CandidateMatrix::new(&haplotype_variants, &haplotypes).unwrap();

        //model computation
        let normalization = false;
        let model = Model::new(Likelihood::new(normalization), Prior::new(), Posterior::new());
        let data = Data::new(candidate_matrix, variant_calls);
        let computed_model = model.compute_from_marginal(&Marginal::new(self.max_haplotypes), &data);
        let mut event_posteriors = computed_model.event_posteriors();

        //add variant query and probabilities to the outout table for each event
        let variant_calls: Vec<AlleleFreqDist> = data.variant_calls
            .iter()
            .map(|(_, afd)| afd.clone())
            .collect();
        let mut event_queries: Vec<BTreeMap<VariantID, (AlleleFreq, LogProb)>> = Vec::new();
        if normalization {
            event_posteriors.for_each(|(fractions, _)| {
                let mut vaf_queries: BTreeMap<VariantID, (AlleleFreq, LogProb)> = BTreeMap::new();
                data.candidate_matrix
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
                            }
                            else {
                                denom -= *fraction;
                            }
                        });
                        if denom > NotNan::new(0.0).unwrap() {
                            vaf_sum /= denom;
                        }
                        vaf_sum = NotNan::new((vaf_sum * NotNan::new(100.0).unwrap()).round()).unwrap()
                            / NotNan::new(100.0).unwrap();
                        if !afd.is_empty() && counter > 0 {
                            let answer = afd.vaf_query(&vaf_sum);
                            vaf_queries.insert(*variant_id, (vaf_sum, answer.unwrap()));
                        } else {
                            ()
                        }
                    });
                event_queries.push(vaf_queries);
            });
        } else {
            event_posteriors.for_each(|(fractions, _)| {
                let mut vaf_queries: BTreeMap<VariantID, (AlleleFreq, LogProb)> = BTreeMap::new();
                data.candidate_matrix
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
                            }
                        });
                        if !afd.is_empty() && counter > 0 {
                            let answer = afd.vaf_query(&vaf_sum);
                            vaf_queries.insert(*variant_id, (vaf_sum, answer.unwrap()));
                        } else {
                            ()
                        }
                    });
                event_queries.push(vaf_queries);
            });
        }

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
        let mut event_posteriors = computed_model.event_posteriors(); //compute a second time because event_posteriors can't be cloned from above.
        let (haplotype_frequencies, best_density) = event_posteriors.next().unwrap();
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
        event_posteriors.zip(event_queries.iter().skip(1)).for_each(
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
        //observations: &mut bcf::Reader,
        haplotype_variants: &mut bcf::Reader,
        filtered_ids: &Vec<VariantID>,
        //max_haplotypes: &usize,
    ) -> Result<Self> {
        let mut variant_records = BTreeMap::new();
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
        Ok(HaplotypeVariants(variant_records))
    }

    fn find_plausible_haplotypes(&self, max_haplotypes: usize) -> Result<Self> {
        //prepare the vectors to be used
        let variants: Vec<VariantID> = self.keys().cloned().collect();
        let (_, haplotypes_gt_c) = self.iter().next().unwrap();
        let haplotype_names: Vec<Haplotype> = haplotypes_gt_c.keys().cloned().collect();
        let len_haplotypes = haplotype_names.len();
        //fill the haplotype vector with variant information
        let mut haplotype_vector = Array2::<f32>::zeros((len_haplotypes, self.len()));
        self.iter()
            .enumerate()
            .for_each(|(variant_index, (_, matrices))| {
                matrices
                    .iter()
                    .enumerate()
                    .for_each(|(haplotype_index, (_, (genotype, _)))| {
                        if *genotype {
                            haplotype_vector[[haplotype_index, variant_index]] = 1.0;
                        }
                    });
            });

        // perform k-means clustering
        let dataset = DatasetBase::from(haplotype_vector);
        let rng = thread_rng(); // Random number generator
        let n_clusters = 5;
        let model = KMeans::params_with_rng(n_clusters, rng)
            .max_n_iterations(200)
            .tolerance(1e-5)
            .fit(&dataset)
            .expect("Error while fitting KMeans to the dataset");
        let dataset = model.predict(dataset);

        //store assignment of each haplotype to corresponding clusters and store vectors belonging to their clusters.
        let mut pseudohaplotypes: BTreeMap<i32, Vec<Haplotype>> = BTreeMap::new();
        pseudohaplotypes.insert(0, vec![]);
        pseudohaplotypes.insert(1, vec![]);
        pseudohaplotypes.insert(2, vec![]);
        pseudohaplotypes.insert(3, vec![]);
        pseudohaplotypes.insert(4, vec![]);

        let mut haplotype_vectors_0 = Vec::new();
        let mut haplotype_vectors_1 = Vec::new();
        let mut haplotype_vectors_2 = Vec::new();
        let mut haplotype_vectors_3 = Vec::new();
        let mut haplotype_vectors_4 = Vec::new();

        for (i, haplotype) in haplotype_names.iter().enumerate() {
            if dataset.targets[i] == 0 {
                let mut existing = pseudohaplotypes.get(&0).unwrap().clone();
                existing.push(haplotype.clone());
                pseudohaplotypes.insert(0, existing.to_vec());
                haplotype_vectors_0.push(dataset.records.slice(s![i, 0..self.len()]));
            } else if dataset.targets[i] == 1 {
                let mut existing = pseudohaplotypes.get(&1).unwrap().clone();
                existing.push(haplotype.clone());
                pseudohaplotypes.insert(1, existing.to_vec());
                haplotype_vectors_1.push(dataset.records.slice(s![i, 0..self.len()]));
            } else if dataset.targets[i] == 2 {
                let mut existing = pseudohaplotypes.get(&2).unwrap().clone();
                existing.push(haplotype.clone());
                pseudohaplotypes.insert(2, existing.to_vec());
                haplotype_vectors_2.push(dataset.records.slice(s![i, 0..self.len()]));
            } else if dataset.targets[i] == 3 {
                let mut existing = pseudohaplotypes.get(&3).unwrap().clone();
                existing.push(haplotype.clone());
                pseudohaplotypes.insert(3, existing.to_vec());
                haplotype_vectors_3.push(dataset.records.slice(s![i, 0..self.len()]));
            } else if dataset.targets[i] == 4 {
                let mut existing = pseudohaplotypes.get(&4).unwrap().clone();
                existing.push(haplotype.clone());
                pseudohaplotypes.insert(4, existing.to_vec());
                haplotype_vectors_4.push(dataset.records.slice(s![i, 0..self.len()]));
            }
        }
        dbg!(&pseudohaplotypes);
        dbg!(&haplotype_vectors_0);
        //collect common variants and nonvariants in clusters
        //5726 -> 5629 for cluster_0
        let mut common_variants_indices_0 = Vec::new();
        let mut common_variants_indices_1 = Vec::new();
        let mut common_variants_indices_2 = Vec::new();
        let mut common_variants_indices_3 = Vec::new();
        let mut common_variants_indices_4 = Vec::new();

        for variant_index in 0..dataset.records.shape()[1] {
            let mut temp_haplotypes = Vec::new();
            for haplotype_array in haplotype_vectors_0.iter() {
                temp_haplotypes.push(haplotype_array[variant_index]);
            }
            let any_element = temp_haplotypes[0];
            if temp_haplotypes.iter().all(|h| h == &any_element) {
                common_variants_indices_0.push(variant_index);
            }
        }

        for variant_index in 0..dataset.records.shape()[1] {
            let mut temp_haplotypes = Vec::new();
            for haplotype_array in haplotype_vectors_1.iter() {
                temp_haplotypes.push(haplotype_array[variant_index]);
            }
            let any_element = temp_haplotypes[0];
            if temp_haplotypes.iter().all(|h| h == &any_element) {
                common_variants_indices_1.push(variant_index);
            }
        }

        for variant_index in 0..dataset.records.shape()[1] {
            let mut temp_haplotypes = Vec::new();
            for haplotype_array in haplotype_vectors_2.iter() {
                temp_haplotypes.push(haplotype_array[variant_index]);
            }
            let any_element = temp_haplotypes[0];
            if temp_haplotypes.iter().all(|h| h == &any_element) {
                common_variants_indices_2.push(variant_index);
            }
        }

        for variant_index in 0..dataset.records.shape()[1] {
            let mut temp_haplotypes = Vec::new();
            for haplotype_array in haplotype_vectors_3.iter() {
                temp_haplotypes.push(haplotype_array[variant_index]);
            }
            let any_element = temp_haplotypes[0];
            if temp_haplotypes.iter().all(|h| h == &any_element) {
                common_variants_indices_3.push(variant_index);
            }
        }

        for variant_index in 0..dataset.records.shape()[1] {
            let mut temp_haplotypes = Vec::new();
            for haplotype_array in haplotype_vectors_4.iter() {
                temp_haplotypes.push(haplotype_array[variant_index]);
            }
            let any_element = temp_haplotypes[0];
            if temp_haplotypes.iter().all(|h| h == &any_element) {
                common_variants_indices_4.push(variant_index);
            }
        }

        let mut final_variants_indices: Vec<usize> = Vec::new();
        final_variants_indices.extend(common_variants_indices_0);
        final_variants_indices.extend(common_variants_indices_1);
        final_variants_indices.extend(common_variants_indices_2);
        final_variants_indices.extend(common_variants_indices_3);
        final_variants_indices.extend(common_variants_indices_4);

        let final_variants_indices: Vec<usize> =
            final_variants_indices.into_iter().unique().collect();
        dbg!(&final_variants_indices.len());
        let mut pseudohaplotypes_variants = BTreeMap::new();
        for index in final_variants_indices {
            let mut matrix_map: BTreeMap<Haplotype, (bool, bool)> = BTreeMap::new();

            for pseudohaplotype in haplotype_vectors_0.iter() {
                //let pseudohaplotype_name = format!("{}{:?}","pseudohaplotype_", cluster_index);
                let pseudohaplotype_name = Haplotype("Pseudohaplotype_0".to_string());
                matrix_map.insert(
                    pseudohaplotype_name,
                    (pseudohaplotype[index] == 1.0, pseudohaplotype[index] == 1.0),
                );
            }
            for pseudohaplotype in haplotype_vectors_1.iter() {
                //let pseudohaplotype_name = format!("{}{:?}","pseudohaplotype_", cluster_index);
                let pseudohaplotype_name = Haplotype("Pseudohaplotype_1".to_string());
                matrix_map.insert(
                    pseudohaplotype_name,
                    (pseudohaplotype[index] == 1.0, pseudohaplotype[index] == 1.0),
                );
            }
            for pseudohaplotype in haplotype_vectors_2.iter() {
                //let pseudohaplotype_name = format!("{}{:?}","pseudohaplotype_", cluster_index);
                let pseudohaplotype_name = Haplotype("Pseudohaplotype_2".to_string());
                matrix_map.insert(
                    pseudohaplotype_name,
                    (pseudohaplotype[index] == 1.0, pseudohaplotype[index] == 1.0),
                );
            }
            for pseudohaplotype in haplotype_vectors_3.iter() {
                //let pseudohaplotype_name = format!("{}{:?}","pseudohaplotype_", cluster_index);
                let pseudohaplotype_name = Haplotype("Pseudohaplotype_3".to_string());
                matrix_map.insert(
                    pseudohaplotype_name,
                    (pseudohaplotype[index] == 1.0, pseudohaplotype[index] == 1.0),
                );
            }
            for pseudohaplotype in haplotype_vectors_4.iter() {
                //let pseudohaplotype_name = format!("{}{:?}","pseudohaplotype_", cluster_index);
                let pseudohaplotype_name = Haplotype("Pseudohaplotype_4".to_string());
                matrix_map.insert(
                    pseudohaplotype_name,
                    (pseudohaplotype[index] == 1.0, pseudohaplotype[index] == 1.0),
                );
            }
            pseudohaplotypes_variants.insert(variants[index].clone(), matrix_map);
        }


        // let model = Model::new(Likelihood::new(), Prior::new(), Posterior::new());
        // let data = Data::new(candidate_matrix, variant_calls);
        // let mut event_posteriors = model.compute_from_marginal(&Marginal::new(n_haplotypes), &data).event_posteriors();
        // let (haplotype_fractions, _) = event_posteriors.next().unwrap();
        
        dbg!(&pseudohaplotypes_variants);
        Ok(HaplotypeVariants(pseudohaplotypes_variants))
    }
}
#[derive(Debug, Clone, Derefable)]
pub(crate) struct AlleleFreqDist(#[deref] BTreeMap<AlleleFreq, f64>);

impl AlleleFreqDist {
    pub(crate) fn vaf_query(&self, vaf: &AlleleFreq) -> Option<LogProb> {
        if self.contains_key(&vaf) {
            Some(LogProb::from(PHREDProb(*self.get(&vaf).unwrap())))
        } else {
            let (x_0, y_0) = self.range(..vaf).next_back().unwrap();
            let (x_1, y_1) = self.range(vaf..).next().unwrap();
            let density =
                NotNan::new(*y_0).unwrap() + (*vaf - *x_0) * (*y_1 - *y_0) / (*x_1 - *x_0); //calculation of density for given vaf by linear interpolation
            Some(LogProb::from(PHREDProb(NotNan::into_inner(density))))
        }
    }
}

#[derive(Derefable, Debug)]
pub(crate) struct CandidateMatrix(#[deref] BTreeMap<VariantID, (BitVec, BitVec)>);

impl CandidateMatrix {
    pub(crate) fn new(
        haplotype_variants: &HaplotypeVariants,
        haplotypes: &Vec<Haplotype>,
    ) -> Result<Self> {
        let mut candidate_matrix = BTreeMap::new();
        haplotype_variants.iter().for_each(|(variant_id, bmap)| {
            let mut haplotype_variants_gt = BitVec::new();
            let mut haplotype_variants_c = BitVec::new();
            bmap.iter().for_each(|(_, (gt, c))| {
                haplotype_variants_gt.push(*gt);
                haplotype_variants_c.push(*c);
            });
            candidate_matrix.insert(*variant_id, (haplotype_variants_gt, haplotype_variants_c));
        });
        Ok(CandidateMatrix(candidate_matrix))
    }
}

#[derive(Derefable, DerefMut, Debug, Clone)]
pub(crate) struct VariantCalls(#[deref] BTreeMap<VariantID, AlleleFreqDist>);

impl VariantCalls {
    pub(crate) fn new(variant_calls: &mut bcf::Reader) -> Result<Self> {
        let mut calls = BTreeMap::new();
        for record_result in variant_calls.records() {
            let mut record = record_result?;
            record.unpack();
            let afd_utf = record.format(b"AFD").string()?;
            let afd = std::str::from_utf8(afd_utf[0]).unwrap();
            let read_depths = record.format(b"DP").integer().unwrap();
            if read_depths[0] != &[0] {
                //because some afd strings are just "." and that throws an error while splitting below.
                let variant_id: i32 = String::from_utf8(record.id())?.parse().unwrap();
                let mut vaf_density = BTreeMap::new();
                for pair in afd.split(',') {
                    if let Some((vaf, density)) = pair.split_once("=") {
                        let (vaf, density): (AlleleFreq, f64) =
                            (vaf.parse().unwrap(), density.parse().unwrap());
                        vaf_density.insert(vaf, density);
                    }
                }
                calls.insert(VariantID(variant_id), AlleleFreqDist(vaf_density));
            }
        }
        Ok(VariantCalls(calls))
    }
}

// fn (num_variants: i32, haplotype_vectors: &Vec<Array>) -> {
//     for variant_index in 0..num_variants {
//         let mut temp_haplotypes = Vec::new();
//         for haplotype_array in haplotype_vectors_0.iter() {
//             temp_haplotypes.push(haplotype_array[variant_index]);
//         }
//         let any_element = temp_haplotypes[0];
//         if temp_haplotypes.iter().all(|h| h == &any_element) {
//             common_variants_0.push(variants[variant_index]);
//         }
//     }
// }
