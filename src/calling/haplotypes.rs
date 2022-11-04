use crate::model::{AlleleFreq, Data, HaplotypeFractions, Likelihood, Marginal, Posterior, Prior};
use anyhow::Result;
use bio::stats::{bayesian::model::Model, probs::LogProb, PHREDProb, Prob};
use bv::BitVec;
use derefable::Derefable;
use derive_builder::Builder;
use derive_deref::DerefMut;
use derive_new::new;
use good_lp::IntoAffineExpression;
use good_lp::*;
use good_lp::{variable, Expression};
use linfa::prelude::*;
use linfa_clustering::KMeans;
use ndarray::prelude::*;
use ordered_float::NotNan;
use rand::prelude::*;
use rand_xoshiro::Xoshiro256Plus;
use rust_htslib::bcf::{self, record::GenotypeAllele::Unphased, Read};
use serde::Serialize;
use serde_json::json;
use std::collections::{BTreeMap, HashMap};
use std::error::Error;
use std::fs::File;
use std::path::Path;
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
        let variant_calls = VariantCalls::new(&mut self.variant_calls)?;
        let variant_ids: Vec<VariantID> = variant_calls.keys().cloned().collect();
        //dbg!(&variant_ids);
        let mut haplotype_variants =
            HaplotypeVariants::new(&mut self.haplotype_variants, &variant_ids)?;
        let mut new_haplotype_variants: BTreeMap<VariantID,BTreeMap<Haplotype, (VariantStatus, bool)>> = BTreeMap::new();
        // if you're only interested in seeing some haplotypes.
        // haplotype_variants.iter().for_each(|(variant, haplotypes)| {
        //     let filtered = haplotypes.iter().filter(|(haplotype, (variant_status, covered))|haplotype.to_string() == "A*03:01".to_string() || haplotype.to_string() == "A*24:02".to_string() || haplotype.to_string() == "A*03:36N".to_string()).map(|(haplotype,(variant_status, covered))|(haplotype.clone(), (variant_status.clone(),*covered))).collect::<BTreeMap<Haplotype, (VariantStatus, bool)>>();
        //     new_haplotype_variants.insert(variant.clone(), filtered);
        // });
        // let new_haplotype_variants = HaplotypeVariants(new_haplotype_variants);

        // let haplotype_variants =
        //     haplotype_variants.find_plausible_haplotypes(&variant_calls, self.max_haplotypes)?;

        let (_, haplotype_matrix) = haplotype_variants.iter().next().unwrap();
        let haplotypes: Vec<Haplotype> = haplotype_matrix.keys().cloned().collect();
        dbg!(&haplotypes);
        let candidate_matrix = CandidateMatrix::new(&haplotype_variants).unwrap();
        linear_program(&candidate_matrix, &haplotypes, &variant_calls);

        //model computation
        // let upper_bond = NotNan::new(1.0).unwrap();
        // let model = Model::new(Likelihood::new(), Prior::new(), Posterior::new());
        // let data = Data::new(candidate_matrix, variant_calls);
        // let computed_model =
        //     model.compute_from_marginal(&Marginal::new(self.max_haplotypes, upper_bond), &data);
        // let event_posteriors = computed_model.event_posteriors();

        // //plotting
        // event_posteriors
        //     .enumerate()
        //     .for_each(|(event_i, (fractions, _))| {
        //         let json = include_str!("../../templates/fractions_barchart.json");
        //         let mut blueprint: serde_json::Value = serde_json::from_str(json).unwrap();
        //         let mut plot_data_variants = Vec::new();
        //         let mut plot_data_haplotype_variants = Vec::new();
        //         let mut plot_data_haplotype_fractions = Vec::new();

        //         let candidate_matrix_values: Vec<(Vec<VariantStatus>, BitVec)> =
        //             data.candidate_matrix.values().cloned().collect();
        //         // let mut final_prob = LogProb::ln_one();
        //         candidate_matrix_values
        //             .iter()
        //             .zip(data.variant_calls.iter())
        //             .for_each(|((genotypes, covered), (variant_id, (af, afd)))| {
        //                 let mut denom = NotNan::new(1.0).unwrap();
        //                 let mut vaf_sum = NotNan::new(0.0).unwrap();
        //                 fractions
        //                     .iter()
        //                     .zip(haplotypes.iter())
        //                     .enumerate()
        //                     .for_each(|(i, (fraction, haplotype))| {
        //                         if genotypes[i] == VariantStatus::Present && covered[i as u64] {
        //                             // vaf_sum += *fraction;
        //                             plot_data_haplotype_fractions.push(
        //                                 dataset_haplotype_fractions {
        //                                     haplotype: haplotype.clone(),
        //                                     fraction: *fraction,
        //                                 },
        //                             );
        //                             plot_data_haplotype_variants.push(dataset_haplotype_variants {
        //                                 variant: *variant_id,
        //                                 haplotype: haplotype.clone(),
        //                             });
        //                             plot_data_variants.push(dataset_variants {
        //                                 variant: *variant_id,
        //                                 vaf: *af,
        //                             });
        //                         } else if genotypes[i] == VariantStatus::Unknown {
        //                             ()
        //                         } else if genotypes[i] == VariantStatus::NotPresent
        //                             && covered[i as u64] == false
        //                         {
        //                             denom -= *fraction;
        //                         }
        //                     });
        //                 if denom > NotNan::new(0.0).unwrap() {
        //                     vaf_sum /= denom;
        //                 }
        //                 //to overcome a bug that results in larger than 1.0 VAF. After around 10 - 15th decimal place, the value becomes larger.
        //                 //In any case, for a direct query to the AFD VAFs (they contain 2 decimal places).
        //                 vaf_sum = NotNan::new((vaf_sum * NotNan::new(100.0).unwrap()).round())
        //                     .unwrap()
        //                     / NotNan::new(100.0).unwrap();
        //             });
        //         //for each event, a plot is generated.
        //         let plot_data_variants = json!(plot_data_variants);
        //         let plot_data_haplotype_variants = json!(plot_data_haplotype_variants);
        //         let plot_data_haplotype_fractions = json!(plot_data_haplotype_fractions);

        //         blueprint["datasets"]["variants"] = plot_data_variants;
        //         blueprint["datasets"]["haplotype_variants"] = plot_data_haplotype_variants;
        //         blueprint["datasets"]["haplotype_fractions"] = plot_data_haplotype_fractions;

        //         let event_id = event_i.to_string();
        //         let output = format!("{}{}{}", "event_", event_id, ".json");
        //         let path = Path::new(&output);
        //         let file = File::create(path).unwrap();
        //         serde_json::to_writer(file, &blueprint);
        //     });
        //plotting

        //add variant query and probabilities to the outout table for each event
        // let variant_calls: Vec<AlleleFreqDist> = data
        //     .variant_calls
        //     .iter()
        //     .map(|(_, (_, afd))| afd.clone())
        //     .collect();
        // let mut event_queries: Vec<BTreeMap<VariantID, (AlleleFreq, LogProb)>> = Vec::new();
        // event_posteriors.for_each(|(fractions, _)| {
        //     let mut vaf_queries: BTreeMap<VariantID, (AlleleFreq, LogProb)> = BTreeMap::new();
        //     data.candidate_matrix
        //         .iter()
        //         .zip(variant_calls.iter())
        //         .for_each(|((variant_id, (genotypes, covered)), afd)| {
        //             let mut denom = NotNan::new(1.0).unwrap();
        //             let mut vaf_sum = NotNan::new(0.0).unwrap();
        //             let mut counter = 0;
        //             fractions.iter().enumerate().for_each(|(i, fraction)| {
        //                 if genotypes[i as u64] == VariantStatus::Present && covered[i as u64] {
        //                     vaf_sum += *fraction;
        //                 } else if genotypes[i as u64] == VariantStatus::Unknown && covered[i as u64]
        //                 {
        //                     todo!();
        //                 } else if genotypes[i as u64] == VariantStatus::Unknown
        //                     && covered[i as u64] == false
        //                 {
        //                     todo!();
        //                 } else if genotypes[i as u64] == VariantStatus::NotPresent
        //                     && covered[i as u64] == false
        //                 {
        //                     denom -= *fraction;
        //                 }
        //             });
        //             if denom > NotNan::new(0.0).unwrap() {
        //                 vaf_sum /= denom;
        //             }
        //             vaf_sum = NotNan::new((vaf_sum * NotNan::new(100.0).unwrap()).round()).unwrap()
        //                 / NotNan::new(100.0).unwrap();
        //             if !afd.is_empty() && counter > 0 {
        //                 let answer = afd.vaf_query(&vaf_sum);
        //                 vaf_queries.insert(*variant_id, (vaf_sum, answer.unwrap()));
        //             } else {
        //                 ()
        //             }
        //         });
        //     event_queries.push(vaf_queries);
        // });

        // // Step 4: print TSV table with results
        // // TODO use csv crate
        // // Columns: posterior_prob, haplotype_a, haplotype_b, haplotype_c, ...
        // // with each column after the first showing the fraction of the respective haplotype
        // let mut wtr = csv::Writer::from_path(self.outcsv.as_ref().unwrap())?;
        // let mut headers: Vec<_> = vec!["density".to_string(), "odds".to_string()];
        // let haplotypes_str: Vec<String> = haplotypes.iter().map(|h| h.to_string()).collect();
        // headers.extend(haplotypes_str);
        // let variant_names = event_queries[0]
        //     .keys()
        //     .map(|key| format!("{:?}", key))
        //     .collect::<Vec<String>>();
        // headers.extend(variant_names); //add variant names as separate columns
        // wtr.write_record(&headers)?;

        // //write best record on top
        // let mut records = Vec::new();
        // let mut event_posteriors = computed_model.event_posteriors(); //compute a second time because event_posteriors can't be cloned from above.
        // let (haplotype_frequencies, best_density) = event_posteriors.next().unwrap();
        // let best_odds = 1;
        // let format_f64 = |number: f64, records: &mut Vec<String>| {
        //     if number <= 0.01 {
        //         records.push(format!("{:+.2e}", number))
        //     } else {
        //         records.push(format!("{:.2}", number))
        //     }
        // };
        // format_f64(best_density.exp(), &mut records);
        // records.push(best_odds.to_string());
        // let format_freqs = |frequency: NotNan<f64>, records: &mut Vec<String>| {
        //     if frequency <= NotNan::new(0.01).unwrap() {
        //         records.push(format!("{:+.2e}", NotNan::into_inner(frequency)))
        //     } else {
        //         records.push(format!("{:.2}", frequency))
        //     }
        // };
        // haplotype_frequencies
        //     .iter()
        //     .for_each(|frequency| format_freqs(*frequency, &mut records));
        // //add vaf queries and probabilities for the first event to the output table
        // let queries: Vec<(AlleleFreq, LogProb)> = event_queries
        //     .iter()
        //     .next()
        //     .unwrap()
        //     .values()
        //     .cloned()
        //     .collect();
        // queries.iter().for_each(|(query, answer)| {
        //     let prob = f64::from(Prob::from(*answer));
        //     if prob <= 0.01 {
        //         records.push(format!("{}{}{:+.2e}", query, ":", prob));
        //     } else {
        //         records.push(format!("{}{}{:.2}", query, ":", prob));
        //     }
        // });
        // wtr.write_record(records)?;

        // //write the rest of the records
        // event_posteriors.zip(event_queries.iter().skip(1)).for_each(
        //     |((haplotype_frequencies, density), queries)| {
        //         let mut records = Vec::new();
        //         let odds = (density - best_density).exp();
        //         format_f64(density.exp(), &mut records);
        //         format_f64(odds, &mut records);
        //         haplotype_frequencies
        //             .iter()
        //             .for_each(|frequency| format_freqs(*frequency, &mut records));

        //         queries.iter().for_each(|(_, (query, answer))| {
        //             let prob = f64::from(Prob::from(*answer));
        //             if prob <= 0.01 {
        //                 records.push(format!("{}{}{:+.2e}", query, ":", prob));
        //             } else {
        //                 records.push(format!("{}{}{:.2}", query, ":", prob));
        //             }
        //         });
        //         wtr.write_record(records).unwrap();
        //     },
        // );
        Ok(())
    }
}

#[derive(Serialize, Debug)]
pub(crate) struct dataset_variants {
    variant: VariantID,
    vaf: f32,
}
#[derive(Serialize, Debug)]
pub(crate) struct dataset_haplotype_variants {
    variant: VariantID,
    haplotype: Haplotype,
}
#[derive(Serialize, Debug)]
pub(crate) struct dataset_haplotype_fractions {
    haplotype: Haplotype,
    fraction: AlleleFreq,
}

#[derive(Derefable, Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize)]
pub(crate) struct Haplotype(#[deref] String);

#[derive(Derefable, Debug, Copy, Clone, PartialEq, Eq, Hash, Ord, PartialOrd, Serialize)]
pub(crate) struct VariantID(#[deref] i32);

#[derive(Clone, Debug, Eq, PartialEq, PartialOrd)]
pub enum VariantStatus {
    Present,
    NotPresent,
    Unknown,
}

#[derive(Derefable, Debug, Clone, PartialEq, Eq, PartialOrd)]
pub(crate) struct HaplotypeVariants(
    #[deref] BTreeMap<VariantID, BTreeMap<Haplotype, (VariantStatus, bool)>>,
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
                        if *gta == Unphased(1) {
                            matrices.insert(
                                haplotype.clone(),
                                (VariantStatus::Present, loci[index] == &[1]),
                            );
                        } else {
                            matrices.insert(
                                haplotype.clone(),
                                (VariantStatus::NotPresent, loci[index] == &[1]),
                            );
                        }
                    }
                }
                variant_records.insert(variant_id, matrices);
            }
        }
        Ok(HaplotypeVariants(variant_records))
    }

    fn find_plausible_haplotypes(
        &self,
        variant_calls: &VariantCalls,
        max_haplotypes: usize,
    ) -> Result<Self> {
        //apply clustering and assign haplotypes to pseudohaplotypes
        let upper_bond = NotNan::new(1.0).unwrap();
        let (variants, pseudohaplotypes, haplotype_fractions) =
            self.cluster_and_run_model(variant_calls, max_haplotypes, &upper_bond)?;

        //TODO: recursively cluster and run the model resulting in the actual combination of haplotypes.
        self.recursive_clustering(
            &variant_calls,
            &haplotype_fractions,
            &pseudohaplotypes,
            max_haplotypes,
        );

        Ok(self.clone())
    }

    fn recursive_clustering(
        &self,
        variant_calls: &VariantCalls,
        fractions: &HaplotypeFractions,
        pseudohaplotypes: &BTreeMap<usize, Vec<Haplotype>>,
        max_haplotypes: usize,
    ) -> Result<()> {
        let mut variants: Vec<VariantID> = Vec::new();
        for (index, fraction) in fractions.iter().enumerate() {
            let haplotypes_in_cluster = pseudohaplotypes.get(&index).unwrap();
            if fraction > &NotNan::new(0.00).unwrap() {
                if haplotypes_in_cluster.len() > max_haplotypes {
                    let selected_haplotypes = haplotypes_in_cluster.clone();
                    //dbg!(&selected_haplotypes);
                    let haplotype_variants_selected =
                        self.filter_haplotypes(&selected_haplotypes).unwrap();
                    //apply clustering and assign haplotypes to pseudohaplotypes
                    let (variants, pseudohaplotypes, haplotype_fractions) =
                        haplotype_variants_selected.cluster_and_run_model(
                            variant_calls,
                            max_haplotypes,
                            fraction,
                        )?;
                    // dbg!(&pseudohaplotypes);
                    // dbg!(&haplotype_fractions);
                    self.recursive_clustering(
                        variant_calls,
                        &haplotype_fractions,
                        &pseudohaplotypes,
                        max_haplotypes,
                    );
                } else {
                    let filtered_haplotype_variants: BTreeMap<
                        VariantID,
                        BTreeMap<Haplotype, (VariantStatus, bool)>,
                    > = self
                        .iter()
                        .filter(|&(k, _)| variants.contains(k))
                        .map(|(k, v)| (k.clone(), v.clone()))
                        .collect();
                    let filtered_haplotype_variants =
                        HaplotypeVariants(filtered_haplotype_variants);
                    let filtered_variant_calls: BTreeMap<VariantID, (f32, AlleleFreqDist)> =
                        variant_calls
                            .iter()
                            .filter(|&(k, (_, _))| variants.contains(k))
                            .map(|(k, v)| (k.clone(), v.clone()))
                            .collect();
                    let filtered_variant_calls = VariantCalls(filtered_variant_calls);

                    let selected_haplotypes = haplotypes_in_cluster.clone();
                    // dbg!(&selected_haplotypes);
                    let haplotype_variants_selected = filtered_haplotype_variants
                        .filter_haplotypes(&selected_haplotypes)
                        .unwrap();

                    let model = Model::new(Likelihood::new(), Prior::new(), Posterior::new());
                    let candidate_matrix =
                        CandidateMatrix::new(&haplotype_variants_selected.clone()).unwrap();
                    let data = Data::new(candidate_matrix, filtered_variant_calls.clone());
                    let computed_model = model.compute_from_marginal(
                        &Marginal::new(selected_haplotypes.len(), *fraction),
                        &data,
                    );
                    let mut event_posteriors = computed_model.event_posteriors();
                    let (haplotype_fractions, _) = event_posteriors.next().unwrap();
                    // dbg!(&haplotypes_in_cluster);
                    // dbg!(&haplotype_fractions);
                }
            }
        }
        Ok(())
    }

    fn cluster_and_run_model(
        &self,
        variant_calls: &VariantCalls,
        max_haplotypes: usize,
        upper_bond: &NotNan<f64>,
    ) -> Result<(
        Vec<VariantID>,
        BTreeMap<usize, Vec<Haplotype>>,
        HaplotypeFractions,
    )> {
        //STEP 1: cluster haplotypes
        //prepare the vectors to be used
        let variants: Vec<VariantID> = self.keys().cloned().collect();
        let (_, haplotypes_gt_c) = self.iter().next().unwrap();
        let haplotype_names: Vec<Haplotype> = haplotypes_gt_c.keys().cloned().collect();
        let len_haplotypes = haplotype_names.len();

        //fill genotype and coverage arrays
        let mut genotype_array = Array2::<f32>::zeros((len_haplotypes, self.len()));
        let mut coverage_array = Array2::<f32>::zeros((len_haplotypes, self.len()));
        self.iter()
            .enumerate()
            .for_each(|(variant_index, (_, matrices))| {
                matrices.iter().enumerate().for_each(
                    |(haplotype_index, (_, (genotype, coverage)))| {
                        if *genotype == VariantStatus::Present {
                            //coverage is automatically true aswell.
                            genotype_array[[haplotype_index, variant_index]] = 1.0;
                            coverage_array[[haplotype_index, variant_index]] = 1.0;
                        } else if *coverage {
                            coverage_array[[haplotype_index, variant_index]] = 1.0;
                        }
                    },
                );
            });

        // perform k-means clustering, seeded for reproducibility
        let seed = 42;
        let rng = Xoshiro256Plus::seed_from_u64(seed);
        let dataset = DatasetBase::from(genotype_array);
        let n_clusters = max_haplotypes;
        let model = KMeans::params_with_rng(n_clusters, rng)
            .max_n_iterations(200)
            .tolerance(1e-5)
            .fit(&dataset)
            .expect("Error while fitting KMeans to the dataset");
        let dataset = model.predict(dataset);

        //initialize vectors
        let mut pseudohaplotypes: BTreeMap<usize, Vec<Haplotype>> = BTreeMap::new();
        for cluster_index in 0..n_clusters {
            pseudohaplotypes.insert(cluster_index, vec![]);
        }

        //collect pseudohaplotype groups
        for (haplotype_index, haplotype) in haplotype_names.iter().enumerate() {
            for cluster_index in 0..n_clusters {
                if dataset.targets[haplotype_index] == cluster_index {
                    let mut existing = pseudohaplotypes.get(&cluster_index).unwrap().clone();
                    existing.push(haplotype.clone());
                    pseudohaplotypes.insert(cluster_index, existing.to_vec());
                }
            }
        }
        // dbg!(&pseudohaplotypes);

        //STEP 2: run the model

        //initialize vectors
        let mut genotype_vectors = Vec::new();
        let mut coverage_vectors = Vec::new();

        for cluster_index in 0..n_clusters {
            genotype_vectors.push(vec![]);
            coverage_vectors.push(vec![]);
        }

        //collect variant vectors of haplotypes to corresponding clusters
        for (haplotype_index, haplotype) in haplotype_names.iter().enumerate() {
            for cluster_index in 0..n_clusters {
                if dataset.targets[haplotype_index] == cluster_index {
                    genotype_vectors[cluster_index]
                        .push(dataset.records.slice(s![haplotype_index, 0..self.len()]));
                    coverage_vectors[cluster_index]
                        .push(coverage_array.slice(s![haplotype_index, 0..self.len()]));
                }
            }
        }

        //create pseudohaplotype names
        let mut pseudohaplotype_names: Vec<Haplotype> = Vec::new();
        for i in 0..n_clusters {
            let name = format!("{}{:?}", "Pseudohaplotype_", i);
            pseudohaplotype_names.push(Haplotype(name));
        }

        //creating final HaplotypeVariants, if one variant is consistenly true or false in one pseudohaplotype group, use Present or NotPresent
        //if not make it UNKNOWN
        let mut pseudohaplotypes_variants = BTreeMap::new();
        for variant_index in 0..variants.len() {
            let mut matrix_map: BTreeMap<Haplotype, (VariantStatus, bool)> = BTreeMap::new();
            for (haplotype_index, (genotype_vector, coverage_vector)) in genotype_vectors
                .iter()
                .zip(coverage_vectors.iter())
                .enumerate()
            {
                let mut coverage_at_index = Vec::new();
                for coverage_array in coverage_vector.iter() {
                    coverage_at_index.push(coverage_array[variant_index]);
                }
                let coverage_info_at_index = coverage_at_index.iter().any(|&x| x == 1.0); //if at least one of them covers, the pseudohaplotype is also said to cover the variant.

                let pseudohaplotype_name = pseudohaplotype_names[haplotype_index].clone();
                let mut variants_at_index = Vec::new();
                for genotype_array in genotype_vector.iter() {
                    variants_at_index.push(genotype_array[variant_index]);
                }
                //now check if variant has the same information across all haplotypes
                //then create HaplotypeVariants accordingly.
                let any_element = variants_at_index[0];
                if variants_at_index.iter().all(|h| h == &any_element) {
                    let random_array = genotype_vector[0];
                    if any_element == 1.0 {
                        matrix_map.insert(
                            pseudohaplotype_name,
                            (VariantStatus::Present, coverage_info_at_index),
                        );
                    } else {
                        matrix_map.insert(
                            pseudohaplotype_name,
                            (VariantStatus::NotPresent, coverage_info_at_index),
                        );
                    }
                } else {
                    matrix_map.insert(
                        pseudohaplotype_name,
                        (VariantStatus::Unknown, coverage_info_at_index),
                    );
                }
            }
            pseudohaplotypes_variants.insert(variants[variant_index].clone(), matrix_map);
        }
        // dbg!(&pseudohaplotypes_variants.len());
        // dbg!(&variant_calls.len());
        //make sure VariantCalls have the same variants
        let variant_calls: BTreeMap<VariantID, (f32, AlleleFreqDist)> = variant_calls
            .iter()
            .filter(|&(k, (_, _))| pseudohaplotypes_variants.contains_key(k))
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        let variant_calls = VariantCalls(variant_calls);

        //model computation, only first round for now
        let model = Model::new(Likelihood::new(), Prior::new(), Posterior::new());
        let candidate_matrix =
            CandidateMatrix::new(&HaplotypeVariants(pseudohaplotypes_variants.clone())).unwrap();
        let data = Data::new(candidate_matrix, variant_calls.clone());
        let computed_model =
            model.compute_from_marginal(&Marginal::new(max_haplotypes, *upper_bond), &data);
        let mut event_posteriors = computed_model.event_posteriors();
        let mut event_posteriors_clone = computed_model.event_posteriors();
        for (fractions, prob) in event_posteriors_clone {
            // dbg!(&fractions);
            // dbg!(&prob.exp());
        }
        let (haplotype_fractions, _) = event_posteriors.next().unwrap();
        let variants: Vec<VariantID> = pseudohaplotypes_variants.keys().cloned().collect();

        Ok((variants, pseudohaplotypes, haplotype_fractions.clone()))
    }

    fn filter_haplotypes(&self, haplotypes: &Vec<Haplotype>) -> Result<Self> {
        //collect first record of self and collect the indices of haplotypes for further filtering in the following.
        let mut haplotype_indices: Vec<usize> = Vec::new();
        if let Some((_, bmap)) = self.iter().next() {
            bmap.iter().enumerate().for_each(|(i, (haplotype, _))| {
                if haplotypes.contains(haplotype) {
                    haplotype_indices.push(i)
                }
            });
        };
        //create filtered matrix with the indices of filtered haplotypes.
        let mut new_variants = BTreeMap::new();
        self.iter().for_each(|(variant_id, bmap)| {
            let mut matrix_map = BTreeMap::new();
            bmap.iter()
                .enumerate()
                .for_each(|(i, (haplotype, matrices))| {
                    haplotype_indices.iter().for_each(|j| {
                        if i == *j {
                            matrix_map.insert(haplotype.clone(), matrices.clone());
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
pub(crate) struct CandidateMatrix(#[deref] BTreeMap<VariantID, (Vec<VariantStatus>, BitVec)>);

impl CandidateMatrix {
    pub(crate) fn new(
        haplotype_variants: &HaplotypeVariants,
        // haplotypes: &Vec<Haplotype>,
    ) -> Result<Self> {
        let mut candidate_matrix = BTreeMap::new();
        haplotype_variants.iter().for_each(|(variant_id, bmap)| {
            let mut haplotype_variants_gt = Vec::new();
            let mut haplotype_variants_c = BitVec::new();
            bmap.iter().for_each(|(haplotype, (gt, c))| {
                haplotype_variants_c.push(*c);
                haplotype_variants_gt.push(gt.clone());
            });
            candidate_matrix.insert(*variant_id, (haplotype_variants_gt, haplotype_variants_c));
        });
        Ok(CandidateMatrix(candidate_matrix))
    }
}

#[derive(Derefable, DerefMut, Debug, Clone)]
pub(crate) struct VariantCalls(#[deref] BTreeMap<VariantID, (f32, AlleleFreqDist)>); //The place of f32 is maximum a posteriori estimate of AF.

impl VariantCalls {
    pub(crate) fn new(variant_calls: &mut bcf::Reader) -> Result<Self> {
        let mut calls = BTreeMap::new();
        for record_result in variant_calls.records() {
            let mut record = record_result?;
            record.unpack();
            let prob_absent = record.info(b"PROB_ABSENT").float().unwrap().unwrap()[0];
            let prob_absent_prob = Prob::from(PHREDProb(prob_absent.into()));
            let afd_utf = record.format(b"AFD").string()?;
            let afd = std::str::from_utf8(afd_utf[0]).unwrap();
            let read_depths = record.format(b"DP").integer().unwrap();
            if read_depths[0] != &[0] && (&prob_absent_prob <= &Prob(0.05) || &prob_absent_prob >= &Prob(0.95)) {
                //because some afd strings are just "." and that throws an error while splitting below.
                let variant_id: i32 = String::from_utf8(record.id())?.parse().unwrap();
                let af = (&*record.format(b"AF").float().unwrap()[0]).to_vec()[0];
                //dbg!(&af);
                let mut vaf_density = BTreeMap::new();
                for pair in afd.split(',') {
                    if let Some((vaf, density)) = pair.split_once("=") {
                        let (vaf, density): (AlleleFreq, f64) =
                            (vaf.parse().unwrap(), density.parse().unwrap());
                        vaf_density.insert(vaf, density);
                    }
                }
                calls.insert(VariantID(variant_id), (af, AlleleFreqDist(vaf_density)));
            }
        }
        Ok(VariantCalls(calls))
    }
}

fn linear_program(
    candidate_matrix: &CandidateMatrix,
    haplotypes: &Vec<Haplotype>,
    variant_calls: &VariantCalls,
) -> Result<(), Box<dyn Error>> {
    //first init the problem
    let mut problem = ProblemVariables::new();
    dbg!(&haplotypes);
    //introduce variables
    let variables: Vec<Variable> =
        problem.add_vector(variable().min(0.0).max(1.0), haplotypes.len());

    //init the constraints
    let mut constraints: Vec<Expression> = Vec::new();

    //execute the objective function to fill up the constraints
    objective_function(
        candidate_matrix,
        haplotypes,
        variant_calls,
        &variables,
        &mut constraints,
    );

    //define temporary variables
    let t_vars: Vec<Variable> = problem.add_vector(variable().min(0.0).max(1.0), constraints.len());

    //create the model to minimise the sum of temporary variables
    let mut sum_tvars = Expression::from_other_affine(0.);
    for t_var in t_vars.iter() {
        sum_tvars += t_var.into_expression();
    }
    let mut model = problem.minimise(sum_tvars.clone()).using(default_solver); // multiple solvers available

    //add a constraint to make sure variables sum up to 1.0.
    let mut sum = Expression::from_other_affine(0.);
    for var in variables.iter() {
        sum += Expression::from_other_affine(var);
    }
    model = model.with(constraint!(sum == 1.0));

    //add the constraints to the model
    for (c, t_var) in constraints.iter().zip(t_vars.iter()) {
        model = model.with(constraint!(t_var >= c.clone()));
        model = model.with(constraint!(t_var >= -c.clone()));
    }
    dbg!(&constraints);

    //make sure again each variable has nonnegative values.
    for variable in variables.iter() {
        model = model.with(constraint!(variable >= 0.0.into_expression()));
    }

    //solve the problem with the default solver, i.e. coin_cbc
    let solution = model.solve().unwrap();

    let mut best_variables = Vec::new();
    //finally, print the variables and the sum
    for (i,(var, haplotype)) in variables.iter().zip(haplotypes.iter()).enumerate() {
        println!("v{}, {}={}",i,haplotype.to_string(),solution.value(*var));
        best_variables.push(solution.value(var.clone()).clone());
    }
    println!("sum = {}", solution.eval(sum_tvars));

    //plot the best result
    dbg!(&best_variables.len());
    plot_solution(&candidate_matrix, &haplotypes, &variant_calls, &best_variables);
    Ok(())
}

fn objective_function(
    candidate_matrix: &CandidateMatrix,
    haplotypes: &Vec<Haplotype>,
    variant_calls: &VariantCalls,
    variables: &Vec<Variable>,
    constraints: &mut Vec<Expression>,
) -> Expression {
    let candidate_matrix_values: Vec<(Vec<VariantStatus>, BitVec)> =
        candidate_matrix.values().cloned().collect();

    //variant-wise iteration
    let mut expr = Expression::from_other_affine(0.); // A constant expression
    for ((genotype_matrix, coverage_matrix), (variant, (af, _))) in
        candidate_matrix_values.iter().zip(variant_calls.iter())
    {
        let mut fraction_cont = Expression::from_other_affine(0.);
        let mut prime_fraction_cont = Expression::from_other_affine(0.);
        let mut vaf = Expression::from_other_affine(0.);
        let mut counter = 0;
        for (i, variable) in variables.iter().enumerate() {
            if coverage_matrix[i as u64] {
                counter += 1;
            }
        }
        if counter == variables.len() {
            dbg!(&variant);
            for (i, (variable, haplotype)) in variables.iter().zip(haplotypes.iter()).enumerate(){
                if genotype_matrix[i] == VariantStatus::Present {
                    fraction_cont += *variable;
                }
            }
            let expr_to_add = fraction_cont - af.clone().into_expression();
            dbg!(&expr_to_add);
            constraints.push(expr_to_add.clone());
            expr += expr_to_add;
        }
    }
    dbg!(&expr);
    expr
}

fn plot_solution(
    candidate_matrix: &CandidateMatrix,
    haplotypes: &Vec<Haplotype>,
    variant_calls: &VariantCalls,
    best_variables: &Vec<f64>,
) -> Result<()> {
    let candidate_matrix_values: Vec<(Vec<VariantStatus>, BitVec)> =
        candidate_matrix.values().cloned().collect();

    //variant-wise iteration
    let mut expr = Expression::from_other_affine(0.); // A constant expression

    //plot
    let json = include_str!("../../templates/fractions_barchart.json");
    let mut blueprint: serde_json::Value = serde_json::from_str(json).unwrap();
    let mut plot_data_variants = Vec::new();
    let mut plot_data_haplotype_variants = Vec::new();
    let mut plot_data_haplotype_fractions = Vec::new();

    for ((genotype_matrix, coverage_matrix), (variant_id, (af, _))) in
        candidate_matrix_values.iter().zip(variant_calls.iter())
    {
        let mut fraction_cont = Expression::from_other_affine(0.);
        let mut prime_fraction_cont = Expression::from_other_affine(0.);
        let mut vaf = Expression::from_other_affine(0.);
        let mut counter = 0;
        for (i, variable) in best_variables.iter().enumerate() {
            if coverage_matrix[i as u64] {
                counter += 1;
            }
        }
        if counter == best_variables.len() {
            for (i, (variable, haplotype)) in best_variables.iter().zip(haplotypes.iter()).enumerate(){
                if genotype_matrix[i] == VariantStatus::Present {
                    fraction_cont += *variable;
                    plot_data_haplotype_fractions.push(
                        dataset_haplotype_fractions {
                            haplotype: haplotype.clone(),
                            fraction: NotNan::new(*variable).unwrap(),
                        },
                    );
                    plot_data_haplotype_variants.push(dataset_haplotype_variants {
                        variant: *variant_id,
                        haplotype: haplotype.clone(),
                    });
                    plot_data_variants.push(dataset_variants {
                        variant: *variant_id,
                        vaf: af.clone(),
                    });
                }
            }
            let expr_to_add = fraction_cont - af.clone().into_expression();
            expr += expr_to_add;
        }
    }

    //for each event, a plot is generated.
    let plot_data_variants = json!(plot_data_variants);
    let plot_data_haplotype_variants = json!(plot_data_haplotype_variants);
    let plot_data_haplotype_fractions = json!(plot_data_haplotype_fractions);

    blueprint["datasets"]["variants"] = plot_data_variants;
    blueprint["datasets"]["haplotype_variants"] = plot_data_haplotype_variants;
    blueprint["datasets"]["haplotype_fractions"] = plot_data_haplotype_fractions;

    let event_id = "best_solution".to_string();
    let output = format!("{}{}{}", "event_", event_id, ".json");
    let path = Path::new(&output);
    let file = File::create(path).unwrap();
    serde_json::to_writer(file, &blueprint);

    Ok(())
}

