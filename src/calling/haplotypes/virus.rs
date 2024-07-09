use crate::calling::haplotypes::haplotypes;
use crate::calling::haplotypes::haplotypes::{
    CandidateMatrix, Haplotype, HaplotypeVariants, PriorTypes, VariantCalls, VariantID,
    VariantStatus,
};
use crate::model::{Data, Likelihood, Marginal, Posterior, Prior};
use anyhow::Result;
use bio::stats::bayesian::model::Model;
use bv::BitVec;
use derive_builder::Builder;
use log::warn;

use ordered_float::NotNan;

use rust_htslib::bcf::{self};

use std::collections::BTreeMap;
use std::fs;
use std::str::FromStr;
use std::{path::PathBuf, str};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    candidates_folder: PathBuf,
    variant_calls: bcf::Reader,
    outcsv: PathBuf,
    prior: String,
    lp_cutoff: f64,
    enable_equivalence_class_constraint: bool,
    threshold_considered_variants: f64,
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

            if variant_calls
                .check_variant_threshold(&haplotype_variants, self.threshold_considered_variants)?
            {
                //filter variants for variants from variant calls
                // let filtered_haplotype_variants:  BTreeMap<VariantID, BTreeMap<Haplotype, (VariantStatus, bool)>> = BTreeMap::new();
                // for (variant,haplotype_map) in haplotype_variants.iter() {
                //     if variant_ids.contains(&variant) {
                //         filtered_haplotype_variants.insert(variant.clone(), haplotype_map.clone());
                //     }
                // }
                // let filtered_haplotype_variants = HaplotypeVariants(filtered_haplotype_variants);

                let filtered_haplotype_variants =
                    haplotype_variants.filter_for_variants(&variant_ids)?;

                let (_, haplotype_matrix) = filtered_haplotype_variants.iter().next().unwrap();
                let haplotypes: Vec<Haplotype> = haplotype_matrix.keys().cloned().collect();
                dbg!(&haplotypes);

                //find the haplotypes to prioritize
                let candidate_matrix = CandidateMatrix::new(&filtered_haplotype_variants).unwrap();

                //currently, ideal variant distance used for extension is 0 for viruses.
                let num_variant_distance: i64 = 1;

                //employ the lineaar program
                let lp_haplotypes = haplotypes::linear_program(
                    &self.outcsv,
                    &candidate_matrix,
                    &haplotypes,
                    &variant_calls,
                    self.lp_cutoff,
                    num_variant_distance,
                )?;
                dbg!(&lp_haplotypes);

                //take only haplotypes that are found by lp
                let lp_haplotype_variants = filtered_haplotype_variants
                    .find_plausible_haplotypes(&variant_calls, &lp_haplotypes)?; //fix: find_plausible haplotypes should only contain the list of "haplotypes" given as parameter
                dbg!(&lp_haplotype_variants);

                //make sure lp_haplotypes sorted the same as in haplotype_variants
                let (_, haplotype_matrix) = lp_haplotype_variants.iter().next().unwrap();
                let final_haplotypes: Vec<Haplotype> = haplotype_matrix.keys().cloned().collect();
                dbg!(&final_haplotypes);

                //construct candidate matrix
                let candidate_matrix = CandidateMatrix::new(&lp_haplotype_variants).unwrap();

                //
                let eq_graph = lp_haplotype_variants
                    .find_equivalence_class("virus")
                    .unwrap();
                dbg!(&eq_graph);

                //1-) model computation for chosen prior
                let prior = PriorTypes::from_str(&self.prior).unwrap();
                let upper_bond = NotNan::new(1.0).unwrap();
                let model = Model::new(
                    Likelihood::new(),
                    Prior::new(prior.clone()),
                    Posterior::new(),
                );
                let data = Data::new(candidate_matrix.clone(), variant_calls.clone());
                let computed_model = model.compute_from_marginal(
                    &Marginal::new(
                        final_haplotypes.len(),
                        final_haplotypes.clone(),
                        upper_bond,
                        prior,
                        eq_graph,
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

                //write to tsv for nonzero densities
                let mut event_posteriors = Vec::new();
                computed_model
                    .event_posteriors()
                    .for_each(|(fractions, logprob)| {
                        if logprob.exp() != 0.0 {
                            event_posteriors.push((fractions.clone(), logprob.clone()));
                        }
                    });
                haplotypes::write_results(
                    &self.outcsv,
                    &data,
                    &event_posteriors,
                    &final_haplotypes,
                    self.prior.clone(),
                    false,
                )?;

                //plot first 10 posteriors of orthanq output
                haplotypes::plot_densities(
                    &self.outcsv,
                    &event_posteriors,
                    &final_haplotypes,
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
