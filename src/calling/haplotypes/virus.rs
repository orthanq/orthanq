use crate::calling::haplotypes::haplotypes;
use crate::calling::haplotypes::haplotypes::get_event_posteriors;
// use crate::calling::haplotypes::haplotypes::HaplotypeGraphVirus;
// use crate::calling::haplotypes::haplotypes::SimilarL;
use crate::calling::haplotypes::haplotypes::{CandidateMatrix, HaplotypeVariants, VariantCalls};

use anyhow::Result;

use derive_builder::Builder;

use ordered_float::NotNan;

use rust_htslib::bcf::{self};

use std::fs;

use std::{path::PathBuf, str};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    haplotype_variants: bcf::Reader,
    variant_calls: bcf::Reader,
    output_folder: PathBuf,
    prior: String,
    lp_cutoff: f64,
    extend_haplotypes: bool,
    num_extend_haplotypes: i64,
    num_constraint_haplotypes: i32,
    output_lp_datavzrd: bool,
    threshold_posterior_density: i32,
}

impl Caller {
    pub fn call(&mut self) -> Result<()> {
        //Step 1: Prepare data and compute the model
        //initially prepare haplotype_variants and variant_calls
        let variant_calls =
            VariantCalls::new(&mut self.variant_calls, &None, &["present".to_string()])?;

        //write blank plots and tsv table if no variants are available.
        if variant_calls.is_empty() {
            //write blank plots, required for the workflow!
            self.output_empty_files()?;
            Ok(())
        } else {
            let haplotype_variants = HaplotypeVariants::new(&mut self.haplotype_variants)?;
            let (event_posteriors, all_haplotypes) = get_event_posteriors(
                &self.output_lp_datavzrd,
                &haplotype_variants,
                &variant_calls,
                "virus",
                &self.prior,
                &self.output_folder,
                self.extend_haplotypes,
                self.num_extend_haplotypes,
                self.num_constraint_haplotypes,
                self.lp_cutoff,
                false,
                None,
                self.threshold_posterior_density,
            )?;

            //find best fractions
            let (best_fractions, _) = event_posteriors.first().unwrap();
            let best_fractions = best_fractions
                .iter()
                .map(|f| NotNan::into_inner(*f))
                .collect::<Vec<f64>>();

            //collect candidate matrix
            let candidate_matrix = CandidateMatrix::new(
                &haplotype_variants
                    .filter_for_haplotypes(&all_haplotypes)
                    .unwrap(),
            )
            .unwrap();
            let candidate_matrix_values: Vec<(bv::BitVec, bv::BitVec)> =
                candidate_matrix.values().cloned().collect();

            //plot best solution
            haplotypes::plot_prediction(
                &self.output_lp_datavzrd,
                &self.output_folder,
                "final",
                &candidate_matrix_values,
                &all_haplotypes,
                &variant_calls,
                &best_fractions,
            )?;

            //write results to tsv
            haplotypes::write_results(
                &self.output_folder.join("predictions.csv"),
                &variant_calls,
                &candidate_matrix,
                &event_posteriors,
                &all_haplotypes,
                true,
            )?;

            //plot first 10 posteriors of orthanq output
            haplotypes::plot_densities(
                &self.output_folder,
                &event_posteriors,
                &all_haplotypes,
                "viral",
                true,
            )?;
            Ok(())
        }
    }
    pub fn output_empty_files(&self) -> Result<()> {
        //write blank plots, required for the workflow!
        fs::create_dir_all(&self.output_folder)?;

        let json: &str = include_str!("../../../templates/final_prediction.json");
        let blueprint: serde_json::Value = serde_json::from_str(json).unwrap();

        for file_name in [
            "lp_solution.json".to_string(),
            "best_solution.json".to_string(),
        ] {
            let blueprint: serde_json::Value = serde_json::from_str(json).unwrap();
            let file = fs::File::create(self.output_folder.join(file_name)).unwrap();
            serde_json::to_writer(file, &blueprint)?;
        }

        //write empty viral solutions
        let file = fs::File::create(self.output_folder.join("viral_solutions.json")).unwrap();
        serde_json::to_writer(file, &blueprint)?;

        //write blank tsv
        let mut wtr = csv::Writer::from_path(self.output_folder.join("predictions.csv"))?;
        let headers: Vec<_> = vec!["density".to_string(), "odds".to_string()];
        wtr.write_record(&headers)?;
        Ok(())
    }
}
