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

use ordered_float::NotNan;

use rust_htslib::bcf::{self};

use std::fs;
use std::str::FromStr;
use std::{path::PathBuf, str};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    haplotype_variants: bcf::Reader,
    variant_calls: bcf::Reader,
    outcsv: PathBuf,
    prior: String,
    lp_cutoff: f64,
}

impl Caller {
    pub fn call(&mut self) -> Result<()> {
        //Step 1: Prepare data and compute the model
        //initially prepare haplotype_variants and variant_calls
        let variant_calls = VariantCalls::new(&mut self.variant_calls)?;

        //write blank plots and tsv table if no variants are available.
        if variant_calls.len() == 0 {
            let mut parent = self.outcsv.clone();
            parent.pop();
            fs::create_dir_all(&parent)?;

            //write blank plots, required for the workflow!
            for file_name in vec![
                "lp_solution.json".to_string(),
                "final_solution.json".to_string(),
            ] {
                let json = include_str!("../../../templates/prediction.json");
                let blueprint: serde_json::Value = serde_json::from_str(json).unwrap();
                let file = fs::File::create(parent.join(file_name)).unwrap();
                serde_json::to_writer(file, &blueprint)?;
            }
            //write blank tsv
            let mut wtr = csv::Writer::from_path(&self.outcsv)?;
            let headers: Vec<_> = vec!["density".to_string(), "odds".to_string()];
            wtr.write_record(&headers)?;
            Ok(())
        } else {
            let variant_ids: Vec<VariantID> = variant_calls.keys().cloned().collect();
            let haplotype_variants =
                HaplotypeVariants::new(&mut self.haplotype_variants, &variant_ids)?;
            let (_, haplotype_matrix) = haplotype_variants.iter().next().unwrap();
            let haplotypes: Vec<Haplotype> = haplotype_matrix.keys().cloned().collect();
            dbg!(&haplotypes);

            //find the haplotypes to prioritize
            let candidate_matrix = CandidateMatrix::new(&haplotype_variants).unwrap();
            let lp_haplotypes = haplotypes::linear_program(
                &self.outcsv,
                &candidate_matrix,
                &haplotypes,
                &variant_calls,
                self.lp_cutoff,
            )?;
            dbg!(&lp_haplotypes);

            //take only haplotypes that are found by lp
            let haplotype_variants =
                haplotype_variants.find_plausible_haplotypes(&variant_calls, &lp_haplotypes)?; //fix: find_plausible haplotypes should only contain the list of "haplotypes" given as parameter
            dbg!(&haplotype_variants);

            //make sure lp_haplotypes sorted the same as in haplotype_variants
            let (_, haplotype_matrix) = haplotype_variants.iter().next().unwrap();
            let final_haplotypes: Vec<Haplotype> = haplotype_matrix.keys().cloned().collect();
            dbg!(&final_haplotypes);

            //construct candidate matrix
            let candidate_matrix = CandidateMatrix::new(&haplotype_variants).unwrap();

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
                &Marginal::new(final_haplotypes.len(), upper_bond, prior),
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

            Ok(())
        }
    }
}
