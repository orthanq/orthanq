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
use std::collections::HashMap;
use std::collections::HashSet;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    haplotype_variants: bcf::Reader,
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
            //prepare haplotype variants
            let haplotype_variants_all = HaplotypeVariants::new(&mut self.haplotype_variants)?;

            //filter variant calls and haplotype variants
            let filtered_calls = variant_calls.without_zero_dp();
            let nonzero_dp_variants: Vec<VariantID> = filtered_calls.keys().cloned().collect();
            let var_filt_haplotype_variants =
                haplotype_variants_all.filter_for_variants(&nonzero_dp_variants)?;

            //find identical haplotypes using variants and prepare LP inputs
            let haplotypes: Vec<Haplotype> = var_filt_haplotype_variants
                .iter()
                .next()
                .unwrap()
                .1
                .keys()
                .cloned()
                .collect();
            let candidate_matrix = CandidateMatrix::new(&var_filt_haplotype_variants).unwrap();
            //generate one map with representative haplotypes as key (required for lp) and one map with all haplotypes as key (required for extension of resulting table)
            let (identical_haplotypes_map_rep, identical_haplotypes_map) =
                candidate_matrix.find_identical_haplotypes(haplotypes);
            let representatives = identical_haplotypes_map_rep.keys().cloned().collect();
            let repr_haplotype_variants =
                var_filt_haplotype_variants.filter_for_haplotypes(&representatives)?;
            let repr_candidate_matrix = CandidateMatrix::new(&repr_haplotype_variants).unwrap();

            //employ the linear program
            //note: extension is disabled at the moment. see notes on linear_program function.
            let lp_haplotypes = haplotypes::linear_program(
                &self.outcsv,
                &repr_candidate_matrix,
                &representatives,
                &filtered_calls,
                self.lp_cutoff,
                // self.extend_haplotypes.unwrap_or(true),
                // self.num_extend_haplotypes, //for now it has to be 0 only
                self.num_constraint_haplotypes,
            )?;

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
                    false, //equivalence class constraint is currently disabled.
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
                        event_posteriors.push((fractions.clone(), logprob.clone()));
                    }
                });

            // Third, extend the table with identical haplotypes
            let (new_event_posteriors, all_haplotypes) = extend_resulting_table(
                &lp_haplotypes,
                &event_posteriors,
                &identical_haplotypes_map,
            )?;

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
fn generate_combinations(
    expanded_fractions: &Vec<AlleleFreq>,
    haplotype_indices: &BTreeMap<Haplotype, usize>,
    identical_haplotypes_map: &BTreeMap<Haplotype, Vec<Haplotype>>,
) -> Vec<Vec<AlleleFreq>> {
    let mut result = vec![expanded_fractions.clone()];

    for (cur_haplotype, &idx_current) in haplotype_indices.iter() {
        let fraction = expanded_fractions[idx_current];

        if fraction > NotNan::new(0.0).unwrap() {
            if let Some(identical_haplotypes) = identical_haplotypes_map.get(cur_haplotype) {
                let mut new_combinations = Vec::new();

                for existing_row in &result {
                    for ident_h in identical_haplotypes {
                        if ident_h != cur_haplotype {
                            let idx_ident_h = haplotype_indices.get(ident_h).unwrap();
                            let mut alt_row = existing_row.clone();
                            alt_row[*idx_ident_h] = fraction;
                            alt_row[idx_current] = NotNan::new(0.0).unwrap();
                            new_combinations.push(alt_row);
                        }
                    }
                }

                result.extend(new_combinations);
            }
        }
    }

    result
}

fn extend_resulting_table(
    representatives: &Vec<Haplotype>,
    event_posteriors: &Vec<(HaplotypeFractions, LogProb)>,
    identical_haplotypes_map: &BTreeMap<Haplotype, Vec<Haplotype>>,
) -> Result<(Vec<(HaplotypeFractions, LogProb)>, Vec<Haplotype>)> {
    let mut all_haplotypes = BTreeSet::new();
    for haplotype in representatives {
        all_haplotypes.insert(haplotype.clone());
        let identical_haplotypes = identical_haplotypes_map.get(haplotype).unwrap();
        all_haplotypes.extend(identical_haplotypes.clone());
    }
    let all_haplotypes: Vec<Haplotype> = all_haplotypes.into_iter().collect();

    let mut new_event_posteriors = Vec::new();

    let haplotype_indices: BTreeMap<Haplotype, usize> = all_haplotypes
        .iter()
        .enumerate()
        .map(|(i, hap)| (hap.clone(), i))
        .collect();

    for (fractions, logprob) in event_posteriors {
        let mut expanded_fractions = vec![NotNan::new(0.0).unwrap(); all_haplotypes.len()];

        for (i, &fraction) in fractions.iter().enumerate() {
            if let Some(&idx) = haplotype_indices.get(&representatives[i]) {
                expanded_fractions[idx] = fraction;
            }
        }

        let all_combinations = generate_combinations(
            &expanded_fractions,
            &haplotype_indices,
            identical_haplotypes_map,
        );

        for combination in all_combinations {
            new_event_posteriors.push((HaplotypeFractions(combination), logprob.clone()));
        }
    }

    Ok((new_event_posteriors, all_haplotypes))
}

fn filter_representatives(
    haplotypes: Vec<Haplotype>,
    distance_matrix: DistanceMatrix,
) -> Vec<Haplotype> {
    let mut representative_set = BTreeSet::new();
    let mut visited = HashSet::new();

    //precompute zero-distance clusters in a HashMap for fast lookup
    let mut zero_distance_clusters: HashMap<Haplotype, Vec<Haplotype>> = HashMap::new();

    for ((h1, h2), &distance) in &*distance_matrix {
        if distance == 0 {
            zero_distance_clusters
                .entry(h1.clone())
                .or_default()
                .push(h2.clone());
            zero_distance_clusters
                .entry(h2.clone())
                .or_default()
                .push(h1.clone());
        }
    }

    for haplotype in &haplotypes {
        if visited.contains(haplotype) {
            continue;
        }

        let mut group = vec![haplotype.clone()];
        visited.insert(haplotype.clone());

        //retrieve precomputed zero-distance neighbors
        if let Some(neighbors) = zero_distance_clusters.get(haplotype) {
            for neighbor in neighbors {
                if !visited.contains(neighbor) {
                    //to prevent redundant addition
                    group.push(neighbor.clone());
                    visited.insert(neighbor.clone());
                }
            }
        }

        //find the smallest haplotype in the group (lexicographically)
        if let Some(min_hap) = group.iter().min() {
            representative_set.insert(min_hap.clone());
        }
    }

    representative_set.into_iter().collect()
}
