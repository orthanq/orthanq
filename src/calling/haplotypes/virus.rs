use crate::calling::haplotypes::haplotypes;
// use crate::calling::haplotypes::haplotypes::find_similar_lp_haplotypes;
// use crate::calling::haplotypes::haplotypes::HaplotypeGraphVirus;
// use crate::calling::haplotypes::haplotypes::SimilarL;
use crate::calling::haplotypes::haplotypes::{
    CandidateMatrix, Haplotype, HaplotypeVariants, PriorTypes, VariantCalls, VariantID,
};

use crate::model::{AlleleFreq, Data, HaplotypeFractions, Likelihood, Marginal, Posterior, Prior};

use anyhow::Result;
use bio::stats::bayesian::model::Model;
use derive_builder::Builder;
use log::warn;

use ordered_float::NotNan;

use rust_htslib::bcf::{self};

use bio::stats::LogProb;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::fs;
use std::str::FromStr;
use std::{path::PathBuf, str};

use super::haplotypes::DistanceMatrix;
use itertools::Itertools;
use std::collections::HashSet;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    candidates_folder: PathBuf,
    variant_calls: bcf::Reader,
    outcsv: PathBuf,
    prior: String,
    lp_cutoff: f64,
    enable_equivalence_class_constraint: bool,
    extend_haplotypes: Option<bool>,
    // threshold_considered_variants: f64,
    threshold_equivalence_class: usize,
    num_extend_haplotypes: i64,
    num_constraint_haplotypes: i32,
}

impl Caller {
    pub fn call(&mut self) -> Result<()> {
        //Step 1: Prepare data and compute the model
        //initially prepare haplotype_variants and variant_calls
        let variant_calls = VariantCalls::new(&mut self.variant_calls)?;

        //write blank plots and tsv table if no variants are available.
        if variant_calls.len() == 0 {
            //write blank plots, required for the workflow!
            self.output_empty_files()?;
            Ok(())
        } else {
            //FIRST, perform linear program using only nonzero DP variants
            //prepare haplotype variants for all variants
            let mut haplotype_variants_rdr =
                bcf::Reader::from_path(self.candidates_folder.join("candidates.vcf"))?;
            let haplotype_variants_all = HaplotypeVariants::new(&mut haplotype_variants_rdr)?;

            //filter variant calls
            let filtered_calls = variant_calls.without_zero_dp();
            let nonzero_dp_variants: Vec<VariantID> = filtered_calls.keys().cloned().collect();
            //filter haplotype variants
            let var_filt_haplotype_variants =
                haplotype_variants_all.filter_for_variants(&nonzero_dp_variants)?;

            //prepare inputs for linear program
            let haplotypes: Vec<Haplotype> = var_filt_haplotype_variants
                .iter()
                .next()
                .unwrap()
                .1
                .keys()
                .cloned()
                .collect();
            let candidate_matrix = CandidateMatrix::new(&var_filt_haplotype_variants).unwrap();
            dbg!(&"ch1");
            //use representative haplotypes for lp
            let distance_matrix_before_lp = var_filt_haplotype_variants
            .find_equivalence_classes_hamming_distance("virus")
            .unwrap();
            let distance_matrix_before_lp = candidate_matrix
            .compute_hamming_distance(&haplotypes);
            dbg!(&"ch2");
            dbg!(&distance_matrix_before_lp);
            let representatives =
            filter_representatives(haplotypes, distance_matrix_before_lp.clone());
            dbg!(&representatives.len());
            dbg!(&representatives);
            let repr_haplotype_variants = var_filt_haplotype_variants.filter_for_haplotypes(&representatives)?;
            let repr_candidate_matrix = CandidateMatrix::new(&repr_haplotype_variants).unwrap();
            //employ the linear program and find the resulting haplotypes that are found and *extended* depending on --extend-haplotypes and --num-extend-haplotypes (0 default)
            let (extended_lp_haplotypes, lp_haplotypes) = haplotypes::linear_program(
                &self.outcsv,
                &repr_candidate_matrix,
                &representatives,
                &filtered_calls,
                self.lp_cutoff,
                self.extend_haplotypes.unwrap_or(true),
                self.num_extend_haplotypes, //for now it has to be 0 only
                self.num_constraint_haplotypes,
            )?;
            dbg!(&lp_haplotypes);
            //SECOND, model evaluation using ALL variants but only the LP- selected haplotypes
            //prepare inputs of model evaluation
            let hap_filt_haplotype_variants =
                haplotype_variants_all.filter_for_haplotypes(&lp_haplotypes)?;
            let model_candidate_matrix = CandidateMatrix::new(&hap_filt_haplotype_variants)?;

            //compute model
            let prior = PriorTypes::from_str(&self.prior).unwrap();
            let upper_bond = NotNan::new(1.0).unwrap();
            let model = Model::new(
                Likelihood::new(),
                Prior::new(prior.clone()),
                Posterior::new(),
            );

            let data = Data::new(model_candidate_matrix.clone(), variant_calls.clone());

            //marginal computation
            let computed_model = model.compute_from_marginal(
                &Marginal::new(
                    lp_haplotypes.len(),
                    lp_haplotypes.clone(),
                    upper_bond,
                    prior,
                    None,
                    self.enable_equivalence_class_constraint,
                    "virus".to_string(),
                ),
                &data,
            );

            //find event posteriors
            let event_posteriors = computed_model.event_posteriors();

            //remove zero densities from the table
            let mut event_posteriors = Vec::new();
            computed_model
                .event_posteriors()
                .for_each(|(fractions, logprob)| {
                    if logprob.exp() != 0.0 {
                        dbg!(&fractions.clone(), &logprob.clone());
                        event_posteriors.push((fractions.clone(), logprob.clone()));
                    }
                });

            //THIRD, extend the table using the distance matrix. Distance matrix contains only the nonzero variants. If all variants are used, there would be no variants that would have 0 distance to each other.
            //extend the resulting table with zero distance haplotypes. For that, compute distance matrix (hamming distance) with lp haplotypes.
            let extended_haplotype_variants =
                var_filt_haplotype_variants.filter_for_haplotypes(&extended_lp_haplotypes)?;

            let distance_matrix = extended_haplotype_variants
                .find_equivalence_classes_hamming_distance("virus")
                .unwrap();

            let (new_event_posteriors, all_haplotypes) =
                extend_resulting_table(&lp_haplotypes, &event_posteriors, &distance_matrix)
                    .unwrap();

            //plot the best solution as final solution plot
            let (best_fractions, _) = new_event_posteriors.iter().next().unwrap();
            let candidate_matrix_all = CandidateMatrix::new(
                &haplotype_variants_all
                    .filter_for_haplotypes(&all_haplotypes)
                    .unwrap(),
            )
            .unwrap()
            .values()
            .cloned()
            .collect();

            let best_fractions = best_fractions
                .iter()
                .map(|f| NotNan::into_inner(*f))
                .collect::<Vec<f64>>();

            haplotypes::plot_prediction(
                &self.outcsv,
                &"final",
                &candidate_matrix_all,
                &all_haplotypes,
                &data.variant_calls,
                &best_fractions,
            )?;

            //write results to tsv
            haplotypes::write_results(
                &self.outcsv,
                &data,
                &new_event_posteriors,
                &all_haplotypes,
                self.prior.clone(),
                false,
            )?;

            //plot first 10 posteriors of orthanq output
            haplotypes::plot_densities(
                &self.outcsv,
                &new_event_posteriors,
                &all_haplotypes,
                "viral",
            )?;
            // } else {
            //     self.output_empty_files()?;
            //     warn!("Insufficient observations from data!");
            // }
            Ok(())
        }
    }
    pub fn output_empty_files(&self) -> Result<()> {
        //write blank plots, required for the workflow!
        let mut parent = self.outcsv.clone();
        parent.pop();
        fs::create_dir_all(&parent)?;

        let json: &str = include_str!("../../../templates/prediction.json");
        let blueprint: serde_json::Value = serde_json::from_str(json).unwrap();

        for file_name in vec![
            "lp_solution.json".to_string(),
            "final_solution.json".to_string(),
        ] {
            let blueprint: serde_json::Value = serde_json::from_str(json).unwrap();
            let file = fs::File::create(parent.join(file_name)).unwrap();
            serde_json::to_writer(file, &blueprint)?;
        }

        //write empty viral solutions
        let file = fs::File::create(parent.join("viral_solutions.json".to_string())).unwrap();
        serde_json::to_writer(file, &blueprint)?;

        //write blank tsv
        let mut wtr = csv::Writer::from_path(&self.outcsv)?;
        let headers: Vec<_> = vec!["density".to_string(), "odds".to_string()];
        wtr.write_record(&headers)?;
        Ok(())
    }
}

// //todo: double check this function
// fn filter_representatives(
//     haplotypes: Vec<Haplotype>,
//     distance_matrix: DistanceMatrix,
// ) -> Vec<Haplotype> {
//     let mut representative_set = BTreeSet::new(); // Use BTreeSet for deterministic ordering
//     let mut visited = HashSet::new();

//     for haplotype in &haplotypes {
//         if visited.contains(haplotype) {
//             continue;
//         }

//         // Collect all haplotypes in the same zero-distance group
//         let mut group = vec![haplotype.clone()];
//         visited.insert(haplotype.clone());

//         for ((h1, h2), &distance) in &*distance_matrix {
//             if distance == 0 {
//                 if h1 == haplotype && !visited.contains(h2) {
//                     group.push(h2.clone());
//                     visited.insert(h2.clone());
//                 } else if h2 == haplotype && !visited.contains(h1) {
//                     group.push(h1.clone());
//                     visited.insert(h1.clone());
//                 }
//             }
//         }

//         // Sort group and pick the smallest haplotype as representative
//         group.sort(); // Ensure consistent ordering
//         representative_set.insert(group[0].clone()); // Lexicographically smallest
//     }

//     representative_set.into_iter().collect() // Convert to Vec while preserving order
// }

fn extend_resulting_table(
    representatives: &Vec<Haplotype>,
    event_posteriors: &Vec<(HaplotypeFractions, LogProb)>,
    distance_matrix: &DistanceMatrix,
) -> Result<(Vec<(HaplotypeFractions, LogProb)>, Vec<Haplotype>)> {
    //initialize new event_posteriors
    let mut new_event_posteriors = Vec::new();

    // collect all lp haplotypes (with all distances)
    let mut all_haplotypes: BTreeSet<Haplotype> = representatives.iter().cloned().collect();
    for ((hap1, hap2), dist) in &**distance_matrix {
        all_haplotypes.insert(hap1.clone());
        all_haplotypes.insert(hap2.clone());
    }
    let all_haplotypes: Vec<Haplotype> = all_haplotypes.into_iter().collect();

    dbg!(&all_haplotypes);

    // create indices for every haplotype
    let haplotype_indices: BTreeMap<Haplotype, usize> = all_haplotypes
        .iter()
        .enumerate()
        .map(|(i, hap)| (hap.clone(), i))
        .collect();
    dbg!(&haplotype_indices);

    for (fractions, logprob) in event_posteriors {
        // create a new fractions vector by extending existing fractions initialized with 0.0
        let mut expanded_fractions: Vec<AlleleFreq> =
            vec![NotNan::new(0.0).unwrap(); all_haplotypes.len()];

        // extend fractions using indices from representative haplotypes
        for (i, &fraction) in fractions.iter().enumerate() {
            if let Some(&idx) = haplotype_indices.get(&representatives[i]) {
                expanded_fractions[idx] = fraction;
            }
        }

        //push the row to event posteriors
        new_event_posteriors.push((
            HaplotypeFractions(expanded_fractions.clone()),
            logprob.clone(),
        ));

        //push alternative solutions to event posteriors
        for (haplotype, idx) in haplotype_indices.iter() {
            // dbg!(&haplotype, &idx);
            let fraction = expanded_fractions[*idx];
            // dbg!(&fraction);
            for ((hap1, hap2), &dist) in &**distance_matrix {
                if dist == 0 && *fraction > 0.0 {
                    // the second condition avoids having duplicate rows
                    if hap1 == haplotype {
                        if let (Some(idx1), Some(idx2)) =
                            (haplotype_indices.get(hap1), haplotype_indices.get(hap2))
                        {
                            let mut alt1 = expanded_fractions.clone();
                            alt1[*idx2] = fraction;
                            alt1[*idx1] = NotNan::new(0.0).unwrap();
                            new_event_posteriors.push((HaplotypeFractions(alt1), logprob.clone()));
                        }
                    } else if hap2 == haplotype {
                        if let (Some(idx1), Some(idx2)) =
                            (haplotype_indices.get(hap1), haplotype_indices.get(hap2))
                        {
                            let mut alt2 = expanded_fractions.clone();
                            alt2[*idx1] = fraction;
                            alt2[*idx2] = NotNan::new(0.0).unwrap();
                            new_event_posteriors.push((HaplotypeFractions(alt2), logprob.clone()));
                        }
                    }
                }
            }
        }
    }
    Ok((new_event_posteriors, all_haplotypes))
}

//todo: double check this function
fn filter_representatives(
    haplotypes: Vec<Haplotype>,
    distance_matrix: DistanceMatrix,
) -> Vec<Haplotype> {
    let mut representative_set = BTreeSet::new(); // Use BTreeSet for deterministic ordering
    let mut visited = HashSet::new();

    for haplotype in &haplotypes {
        if visited.contains(haplotype) {
            continue;
        }

        // Collect all haplotypes in the same zero-distance group
        let mut group = vec![haplotype.clone()];
        visited.insert(haplotype.clone());

        for ((h1, h2), &distance) in &*distance_matrix {
            if distance == 0 {
                if h1 == haplotype && !visited.contains(h2) {
                    group.push(h2.clone());
                    visited.insert(h2.clone());
                } else if h2 == haplotype && !visited.contains(h1) {
                    group.push(h1.clone());
                    visited.insert(h1.clone());
                }
            }
        }

        // Sort group and pick the smallest haplotype as representative
        group.sort(); // Ensure consistent ordering
        representative_set.insert(group[0].clone()); // Lexicographically smallest
    }

    representative_set.into_iter().collect() // Convert to Vec while preserving order
}