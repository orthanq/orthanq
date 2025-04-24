use crate::calling::haplotypes::haplotypes;
use crate::calling::haplotypes::haplotypes::get_event_posteriors;
// use crate::calling::haplotypes::haplotypes::HaplotypeGraphVirus;
// use crate::calling::haplotypes::haplotypes::SimilarL;
use crate::calling::haplotypes::haplotypes::{
    CandidateMatrix, Haplotype, HaplotypeVariants, VariantCalls,
};

use anyhow::Result;

use derive_builder::Builder;

use ordered_float::NotNan;

use rust_htslib::bcf::{self};

use std::collections::BTreeSet;
use std::fs;

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
    extend_haplotypes: bool,
    num_extend_haplotypes: i64,
    num_constraint_haplotypes: i32,
    output_lp_datavzrd: bool
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
            let haplotype_variants = HaplotypeVariants::new(&mut self.haplotype_variants)?;
            let (event_posteriors, all_haplotypes, data) = get_event_posteriors(
                &self.output_lp_datavzrd,

                &haplotype_variants,
                variant_calls,
                &"virus",
                &self.prior,
                &self.outcsv,
                self.extend_haplotypes,
                self.num_extend_haplotypes,
                self.num_constraint_haplotypes,
                self.lp_cutoff,
                false,
                None,
            )?;

            //plot the best solution as final solution plot
            let (best_fractions, _) = event_posteriors.iter().next().unwrap();
            let candidate_matrix = CandidateMatrix::new(
                &haplotype_variants
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
                &self.output_lp_datavzrd,
                &self.outcsv,
                &"final",
                &candidate_matrix,
                &all_haplotypes,
                &data.variant_calls,
                &best_fractions,
            )?;

            //write results to tsv
            haplotypes::write_results(
                &self.outcsv,
                &data,
                &event_posteriors,
                &all_haplotypes,
                self.prior.clone(),
                false,
            )?;

            //plot first 10 posteriors of orthanq output
            haplotypes::plot_densities(&self.outcsv, &event_posteriors, &all_haplotypes, "viral")?;
            Ok(())
        }
    }
    pub fn output_empty_files(&self) -> Result<()> {
        //write blank plots, required for the workflow!
        let mut parent = self.outcsv.clone();
        parent.pop();
        fs::create_dir_all(&parent)?;

        let json: &str = include_str!("../../../templates/final_prediction.json");
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
