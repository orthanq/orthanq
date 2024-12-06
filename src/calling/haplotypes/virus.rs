use crate::calling::haplotypes::haplotypes;
use crate::calling::haplotypes::haplotypes::{
    CandidateMatrix, Haplotype, HaplotypeVariants, PriorTypes, VariantCalls, VariantID,
    VariantStatus,
};
use crate::model::{Data, Likelihood, Marginal, Posterior, Prior, HaplotypeFractions, AlleleFreq};

use anyhow::Result;
use bio::stats::bayesian::model::Model;
use bv::BitVec;
use derive_builder::Builder;
use log::warn;

use ordered_float::NotNan;

use rust_htslib::bcf::{self};

use std::fs;
use std::str::FromStr;
use std::{path::PathBuf, str};
use std::collections::{HashMap,BTreeMap};
use std::collections::{BTreeSet,HashSet};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    candidates_folder: PathBuf,
    variant_calls: bcf::Reader,
    outcsv: PathBuf,
    prior: String,
    lp_cutoff: f64,
    enable_equivalence_class_constraint: bool,
    extend_haplotypes: bool,
    threshold_considered_variants: f64,
    threshold_equivalence_class: usize,
    num_extend_haplotypes: i64,
}

impl Caller {
    pub fn call(&mut self) -> Result<()> {
        //Step 1: Prepare data and compute the model
        //initially prepare haplotype_variants and variant_calls
        let variant_calls = VariantCalls::new(&mut self.variant_calls)?;

        //read candidates vcf
        let haplotype_variants_dir = self.candidates_folder.join("candidates.vcf");
        let mut haplotype_variants_rdr = bcf::Reader::from_path(haplotype_variants_dir)?;

        //write blank plots and tsv table if no variants are available.
        if variant_calls.len() == 0 {

            //write blank plots, required for the workflow!
            self.output_empty_files()?;
            Ok(())
        } else {
            let variant_ids: Vec<VariantID> = variant_calls.keys().cloned().collect();
            let haplotype_variants = HaplotypeVariants::new(&mut haplotype_variants_rdr)?;

            // check if there is enough observations in the data, and do that by checking rate of evaluated variants
            // if this does not pass, print insufficient data 
            if variant_calls
                .check_variant_threshold(&haplotype_variants, self.threshold_considered_variants)?
            {
                // filter haplotype variants to contain only the evaluated variants in the model.
                let filtered_haplotype_variants =
                    haplotype_variants.filter_for_variants(&variant_ids)?;
                
                // output haplotype list and candidate matrix to be used in lp
                let (_, haplotype_matrix) = filtered_haplotype_variants.iter().next().unwrap();
                let haplotypes: Vec<Haplotype> = haplotype_matrix.keys().cloned().collect();
                let candidate_matrix = CandidateMatrix::new(&filtered_haplotype_variants).unwrap();

                // employ the linear program and find the resulting haplotypes that are found and *extended* depending on --extend-haplotypes and --num-extend-haplotypes (0 default)
                let lp_haplotypes = haplotypes::linear_program(
                    &self.outcsv,
                    &candidate_matrix,
                    &haplotypes,
                    &variant_calls,
                    self.lp_cutoff,
                    self.extend_haplotypes,
                    self.num_extend_haplotypes,
                )?;

                //take only haplotypes that are found by lp
                let lp_haplotype_variants = filtered_haplotype_variants
                    .filter_for_haplotypes(&lp_haplotypes)?;

                //compute distance matrix (hamming distance) with lp haplotypes (extended or not)
                let distance_matrix = lp_haplotype_variants
                    .find_equivalence_classes_hamming_distance(
                        "virus",
                    )
                    .unwrap();

                //create a representative haplotype list; haplotypes that come from lp (extended or not) contain haplotypes that are identical by variant set (zero distance - due to insufficient observation in data)
                let representatives = filter_representatives(lp_haplotypes, distance_matrix.clone());
                dbg!(&representatives);

                //then create haplotype variants with only the representative haplotypes
                let final_haplotype_variants = lp_haplotype_variants
                    .filter_for_haplotypes(&representatives)?; 

                //construct candidate matrix
                let candidate_matrix = CandidateMatrix::new(&final_haplotype_variants).unwrap();

                //model computation with representative haplotypes
                let prior = PriorTypes::from_str(&self.prior).unwrap();
                let upper_bond = NotNan::new(1.0).unwrap();
                let model = Model::new(
                    Likelihood::new(),
                    Prior::new(prior.clone()),
                    Posterior::new(),
                );
                let data = Data::new(candidate_matrix.clone(), variant_calls.clone());

                //create a new distance matrix with no entries of 0 distances and use it in marginal computation
                let distance_matrix_nonzero: BTreeMap<(Haplotype, Haplotype), usize> = distance_matrix
                    .iter()
                    .filter(|&(_, &v)| v != 0)
                    .map(|(k, &v)| (k.clone(), v)) 
                    .collect();

                //determine final haplotypes to be used in the model
                let final_haplotypes = representatives.clone();
                dbg!(&final_haplotypes);

                let computed_model = model.compute_from_marginal(
                    &Marginal::new(
                        final_haplotypes.len(),
                        final_haplotypes.clone(),
                        upper_bond,
                        prior,
                        None,
                        Some(distance_matrix_nonzero),
                        self.enable_equivalence_class_constraint,
                        "virus".to_string(),
                    ),
                    &data,
                );
                let mut event_posteriors = computed_model.event_posteriors();
                let (best_fractions, _) = event_posteriors.next().unwrap();

                //Step 2: plot the final solution
                let candidate_matrix_values: Vec<(Vec<VariantStatus>, BitVec)> =
                    data.candidate_matrix.values().cloned().collect();
                let best_fractions = best_fractions
                    .iter()
                    .map(|f| NotNan::into_inner(*f))
                    .collect::<Vec<f64>>();

                haplotypes::plot_prediction(
                    &self.outcsv,
                    &"final",
                    &candidate_matrix_values,
                    &final_haplotypes,
                    &data.variant_calls,
                    &best_fractions,
                )?;

                //remove zero densities from the table
                let mut event_posteriors = Vec::new();
                computed_model
                    .event_posteriors()
                    .for_each(|(fractions, logprob)| {
                        if logprob.exp() != 0.0 {
                            event_posteriors.push((fractions.clone(), logprob.clone()));
                        }
                    });
                dbg!(&final_haplotypes);
                dbg!(&event_posteriors[0]);

                //todo: double check
                // collect all lp haplotypes (with all distances)
                let mut all_haplotypes: HashSet<Haplotype> = representatives.iter().cloned().collect();
                for ((hap1, hap2), dist) in &distance_matrix {
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

                let mut new_event_posteriors = Vec::new();

                for (fractions, logprob) in &event_posteriors {
                    // create a new fractions vector by extending existing fractions initialized with 0.0
                    let mut expanded_fractions: Vec<AlleleFreq> = vec![NotNan::new(0.0).unwrap(); all_haplotypes.len()];

                    // extend fractions using indices from representative haplotypes
                    for (i, &fraction) in fractions.iter().enumerate() {
                        if let Some(&idx) = haplotype_indices.get(&representatives[i]) {
                            expanded_fractions[idx] = fraction;
                        }
                    }

                    //push the row to event posteriors
                    new_event_posteriors.push((HaplotypeFractions(expanded_fractions.clone()), logprob.clone()));

                    //push alternative solutions to event posteriors
                    for (haplotype, idx) in haplotype_indices.iter() {
                        dbg!(&haplotype, &idx);
                        let fraction = expanded_fractions[*idx];
                        dbg!(&fraction);
                        for ((hap1, hap2), &dist) in &distance_matrix {
                            if dist == 0 && fraction > NotNan::new(0.0).unwrap(){ // the second condition avoids having duplicate rows
                                if hap1 == haplotype {
                                    if let (Some(idx1),Some(idx2)) = (haplotype_indices.get(hap1),haplotype_indices.get(hap2)) {
                                        let mut alt1 = expanded_fractions.clone();
                                        alt1[*idx2] = fraction;
                                        alt1[*idx1] = NotNan::new(0.0).unwrap();
                                        new_event_posteriors.push((HaplotypeFractions(alt1), logprob.clone()));
                                    }
                                } else if hap2 == haplotype {
                                    if let (Some(idx1),Some(idx2)) = (haplotype_indices.get(hap1),haplotype_indices.get(hap2)) {
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
                dbg!(&new_event_posteriors.len());
                dbg!(&new_event_posteriors);

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
            } else {
                self.output_empty_files()?;
                warn!("Insufficient observations from data!");
            }
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

//todo: double check this function
fn filter_representatives(
    haplotypes: Vec<Haplotype>,
    distance_matrix: BTreeMap<(Haplotype, Haplotype), usize>,
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

        for ((h1, h2), &distance) in &distance_matrix {
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