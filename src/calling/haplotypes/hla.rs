use crate::calling::haplotypes::haplotypes::{
    AlleleFreqDist, CandidateMatrix, Haplotype, HaplotypeVariants, PriorTypes, VariantCalls,
    VariantID, VariantStatus,
};
use crate::model::{AlleleFreq, Data, HaplotypeFractions, Likelihood, Marginal, Posterior, Prior};
use anyhow::Result;
use bio::stats::{bayesian::model::Model, probs::LogProb, PHREDProb, Prob};
use bv::BitVec;
use core::cmp::Ordering;
use derefable::Derefable;
use derive_builder::Builder;
use derive_deref::DerefMut;

use good_lp::IntoAffineExpression;
use good_lp::*;
use good_lp::{variable, Expression};

use ordered_float::NotNan;
use ordered_float::OrderedFloat;
use quick_xml::events::Event;
use quick_xml::reader::Reader as xml_reader;

use rust_htslib::bcf::{
    self,
    record::GenotypeAllele::{Phased, Unphased},
    Read,
};
use serde::Serialize;
use serde_json::json;
use std::collections::{BTreeMap, HashMap};

use std::fs;
use std::str::FromStr;
use std::{path::PathBuf, str};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    haplotype_variants: bcf::Reader,
    variant_calls: bcf::Reader,
    xml: PathBuf,
    // max_haplotypes: i64,
    // min_norm_counts: f64,
    outcsv: Option<PathBuf>,
    prior: String,
    common_variants: bool,
    lp_cutoff: f64,
}

impl Caller {
    pub fn call(&mut self) -> Result<()> {
        //Step 1: Prepare data and compute the model
        //initially prepare haplotype_variants and variant_calls
        let variant_calls = VariantCalls::new(&mut self.variant_calls)?;

        //write blank plots and tsv table if no variants are available.
        if variant_calls.len() == 0 {
            let mut parent = self.outcsv.clone().unwrap();
            parent.pop();
            fs::create_dir_all(&parent)?;

            //write blank plots, required for the workflow!
            for file_name in vec![
                "lp_solution.json".to_string(),
                "final_solution.json".to_string(),
                "2_field_solutions.json".to_string(),
                "3_field_solutions.json".to_string(),
            ] {
                let json = include_str!("../../../templates/prediction.json");
                let blueprint: serde_json::Value = serde_json::from_str(json).unwrap();
                let file = fs::File::create(parent.join(file_name)).unwrap();
                serde_json::to_writer(file, &blueprint)?;
            }
            //write blank tsv
            let mut wtr = csv::Writer::from_path(self.outcsv.as_ref().unwrap().clone())?;
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

            //check if common_variants is true, if yes use only common variants both for lp and model
            if self.common_variants {
                //include only common variants if common_variants is true
                //this happens by first finding the common variants and then filtering
                //haplotype_variants and variant_calls to only contain those variants
                let common_variants =
                    haplotype_variants.find_common_variants(&variant_calls, &haplotypes)?;
                let haplotype_variants =
                    haplotype_variants.filter_haplotype_variants(&common_variants)?;
                let variant_calls = variant_calls.filter_variant_calls(&common_variants)?;
                dbg!(&common_variants);
                dbg!(&haplotype_variants);
                dbg!(&variant_calls);
            }
            //find the haplotypes to prioritize
            let candidate_matrix = CandidateMatrix::new(&haplotype_variants).unwrap();
            let lp_haplotypes = self.linear_program(
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

            //for diploid-subclonal priors, max N number of haplotypes should be selected
            // let final_haplotypes = final_haplotypes[0..5].to_vec();
            //

            // //kallisto experimentation
            // //prepare KallistoEstimates for only the haplotypes come from LP
            // let kallisto_estimates = KallistoEstimates::new(&self.hdf5_reader, self.min_norm_counts, &final_haplotypes)?;
            // dbg!(&kallisto_estimates);
            // let kallisto_estimates = kallisto_estimates
            //  .select_haplotypes(4)
            //  .unwrap();
            // dbg!(&kallisto_estimates);
            // let kallisto_haplotypes = kallisto_estimates.keys().cloned().collect();
            // let haplotype_variants =
            // haplotype_variants.find_plausible_haplotypes(&variant_calls, &kallisto_haplotypes)?;
            // let (_, haplotype_matrix) = haplotype_variants.iter().next().unwrap();
            // let final_haplotypes: Vec<Haplotype> = haplotype_matrix.keys().cloned().collect();
            // dbg!(&final_haplotypes); //the final ranking of haplotypes
            // let candidate_matrix = CandidateMatrix::new(&haplotype_variants).unwrap();

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
            self.plot_prediction(
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
            //first: 3-field
            self.write_results(
                self.outcsv.as_ref().unwrap().clone(),
                &data,
                &event_posteriors,
                &final_haplotypes,
                self.prior.clone(),
                false,
            )?;
            //second: 2-field
            let (two_field_haplotypes, two_field_event_posteriors) =
                convert_to_two_field(&event_posteriors, &final_haplotypes)?;
            let mut path_for_two_fields =
                PathBuf::from(self.outcsv.as_ref().unwrap().parent().unwrap());
            path_for_two_fields.push("2-field.csv");
            self.write_results(
                path_for_two_fields,
                &data,
                &two_field_event_posteriors,
                &two_field_haplotypes,
                self.prior.clone(),
                false,
            )?;

            //plot first 10 posteriors of orthanq output
            self.plot_densities(&event_posteriors, &final_haplotypes, "3_field")?;
            self.plot_densities(
                &two_field_event_posteriors,
                &two_field_haplotypes,
                "2_field",
            )?;

            //second: convert to G groups
            let mut converted_name = PathBuf::from(self.outcsv.as_ref().unwrap().parent().unwrap());
            converted_name.push("G_groups.csv");
            let allele_to_g_groups = self.convert_to_g().unwrap();
            let mut final_haplotypes_converted: Vec<Haplotype> = Vec::new();
            final_haplotypes.iter().for_each(|haplotype| {
                let mut conv_haplotype = Vec::new();
                allele_to_g_groups.iter().for_each(|(allele, g_group)| {
                    if allele.starts_with(&haplotype.to_string()) {
                        conv_haplotype.push(g_group.to_string());
                    }
                });
                if conv_haplotype.is_empty() {
                    conv_haplotype.push(haplotype.to_string());
                }
                let conv_haplotype = Haplotype(conv_haplotype[0].clone());
                final_haplotypes_converted.push(conv_haplotype);
            });

            self.write_results(
                converted_name,
                &data,
                &event_posteriors,
                &final_haplotypes_converted,
                self.prior.clone(),
                true,
            )?;
            Ok(())
        }
    }
    fn write_results(
        &self,
        out: PathBuf,
        data: &Data,
        event_posteriors: &Vec<(HaplotypeFractions, LogProb)>,
        final_haplotypes: &Vec<Haplotype>,
        _prior: String,
        variant_info: bool,
    ) -> Result<()> {
        //firstly add variant query and probabilities to the outout table for each event
        let variant_calls: Vec<AlleleFreqDist> = data
            .variant_calls
            .iter()
            .map(|(_, (_, afd))| afd.clone())
            .collect();
        let mut event_queries: Vec<BTreeMap<VariantID, (AlleleFreq, LogProb)>> = Vec::new();
        // let event_posteriors = computed_model.event_posteriors();
        if variant_info {
            event_posteriors.iter().for_each(|(fractions, _)| {
                let mut vaf_queries: BTreeMap<VariantID, (AlleleFreq, LogProb)> = BTreeMap::new();
                data.candidate_matrix
                    .iter()
                    .zip(variant_calls.iter())
                    .for_each(|((variant_id, (genotypes, covered)), afd)| {
                        let mut denom = NotNan::new(1.0).unwrap();
                        let mut vaf_sum = NotNan::new(0.0).unwrap();
                        let mut counter = 0;
                        fractions.iter().enumerate().for_each(|(i, fraction)| {
                            if genotypes[i] == VariantStatus::Present && covered[i as u64] {
                                vaf_sum += *fraction;
                                counter += 1;
                            } else if genotypes[i] == VariantStatus::Unknown && covered[i as u64] {
                                todo!();
                            } else if genotypes[i] == VariantStatus::Unknown
                                && covered[i as u64] == false
                            {
                                todo!();
                            } else if genotypes[i] == VariantStatus::NotPresent
                                && covered[i as u64] == false
                            {
                                denom -= *fraction;
                            }
                        });
                        if denom > NotNan::new(0.0).unwrap() {
                            vaf_sum /= denom;
                        }
                        vaf_sum = NotNan::new((vaf_sum * NotNan::new(100.0).unwrap()).round())
                            .unwrap()
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
        }
        // Then,print TSV table with results
        // Columns: posterior_prob, haplotype_a, haplotype_b, haplotype_c, ...
        // with each column after the first showing the fraction of the respective haplotype
        let mut wtr = csv::Writer::from_path(out)?;
        let mut headers: Vec<_> = vec!["density".to_string(), "odds".to_string()];
        let haplotypes_str: Vec<String> = final_haplotypes
            .clone()
            .iter()
            .map(|h| h.to_string())
            .collect();
        headers.extend(haplotypes_str);
        if variant_info {
            let variant_names = event_queries[0]
                .keys()
                .map(|key| format!("{:?}", key))
                .collect::<Vec<String>>();
            headers.extend(variant_names); //add variant names as separate columns
        }
        wtr.write_record(&headers)?;

        //write best record on top
        let mut records = Vec::new();
        // let mut event_posteriors = computed_model.event_posteriors(); //compute a second time because event_posteriors can't be cloned from above.
        let (haplotype_frequencies, best_density) = event_posteriors.iter().next().unwrap();
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

        if variant_info {
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
        }

        wtr.write_record(records)?;

        //write the rest of the records
        dbg!(&event_posteriors);

        if variant_info {
            event_posteriors
                .iter()
                .skip(1)
                .zip(event_queries.iter().skip(1))
                .for_each(|((haplotype_frequencies, density), queries)| {
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
                });
        } else {
            event_posteriors
                .iter()
                .skip(1)
                .for_each(|(haplotype_frequencies, density)| {
                    let mut records = Vec::new();
                    let odds = (density - best_density).exp();
                    format_f64(density.exp(), &mut records);
                    format_f64(odds, &mut records);
                    haplotype_frequencies
                        .iter()
                        .for_each(|frequency| format_freqs(*frequency, &mut records));

                    wtr.write_record(records).unwrap();
                });
        }

        Ok(())
    }
    fn linear_program(
        &self,
        candidate_matrix: &CandidateMatrix,
        haplotypes: &Vec<Haplotype>,
        variant_calls: &VariantCalls,
        lp_cutoff: f64,
    ) -> Result<Vec<Haplotype>> {
        //first init the problem
        let mut problem = ProblemVariables::new();
        //introduce variables
        let variables: Vec<Variable> =
            problem.add_vector(variable().min(0.0).max(1.0), haplotypes.len());

        //init the constraints
        let mut constraints: Vec<Expression> = Vec::new();

        //execute the following function to fill up the constraints and create a haplotype_dict
        let haplotype_dict = collect_constraints_and_variants(
            candidate_matrix,
            haplotypes,
            variant_calls,
            &variables,
            &mut constraints,
        )
        .unwrap();

        //define temporary variables
        let t_vars: Vec<Variable> =
            problem.add_vector(variable().min(0.0).max(1.0), constraints.len());

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
        // dbg!(&constraints);

        //solve the problem with the default solver, i.e. coin_cbc
        let solution = model.solve().unwrap();

        let mut best_variables = Vec::new();
        //finally, print the variables and the sum
        let mut lp_haplotypes = BTreeMap::new();
        for (i, (var, haplotype)) in variables.iter().zip(haplotypes.iter()).enumerate() {
            println!("v{}, {}={}", i, haplotype.to_string(), solution.value(*var));
            best_variables.push(solution.value(var.clone()).clone());
            if solution.value(*var) >= lp_cutoff {
                //the speed of fraction exploration is managable in case of diploid priors
                lp_haplotypes.insert(haplotype.clone(), solution.value(*var).clone());
            }
        }
        println!("sum = {}", solution.eval(sum_tvars));

        //plot the best result
        // dbg!(&best_variables.len());
        let candidate_matrix_values: Vec<(Vec<VariantStatus>, BitVec)> =
            candidate_matrix.values().cloned().collect();
        self.plot_prediction(
            &"lp",
            &candidate_matrix_values,
            &haplotypes,
            &variant_calls,
            &best_variables,
        )?;

        //extend haplotypes found by linear program, add haplotypes that have the same variants to the final list
        //and optionally, sort by hamming distance, take the closest x additional alleles according to 'permitted'
        let mut extended_haplotypes = Vec::new();
        lp_haplotypes.iter().for_each(|(f_haplotype, _)| {
            let variants = haplotype_dict.get(&f_haplotype).unwrap().clone();
            haplotype_dict
                .iter()
                .for_each(|(haplotype, haplotype_variants)| {
                    if &variants == haplotype_variants && !extended_haplotypes.contains(haplotype) {
                        //fix: the last operand '&&' is required to avoid duplicate additions
                        extended_haplotypes.push(haplotype.clone());
                    } else {
                        let permitted: i64 = 3; //this number should be discussed.
                        let mut difference = vec![];
                        for i in haplotype_variants.iter() {
                            if !variants.contains(&i) {
                                difference.push(i);
                            }
                        }
                        if (difference.len() as i64 <= permitted)
                            && ((variants.len() as i64 - haplotype_variants.len() as i64).abs()
                                <= permitted)
                            && !extended_haplotypes.contains(&haplotype)
                        //fix: the last operand '&&' is required to avoid duplicate additions
                        {
                            extended_haplotypes.push(haplotype.clone());
                        }
                    }
                });
        });
        dbg!(&lp_haplotypes);
        //diploid-subclonal max N haplotypes
        // let max_haplotypes = 5;
        let lp_keys: Vec<_> = lp_haplotypes.keys().cloned().collect();
        let lp_values_original: Vec<_> = lp_haplotypes.values().cloned().collect();
        let mut lp_values: Vec<_> = lp_haplotypes.values().cloned().collect();
        lp_values.sort_by(|a, b| OrderedFloat(*b).cmp(&OrderedFloat(*a)));
        let mut selected_haplotypes = Vec::new();
        for x in lp_values.iter() {
            for (i, y) in lp_values_original.iter().enumerate() {
                if x == y {
                    selected_haplotypes.push(lp_keys[i].clone());
                }
            }
        }
        // if max_haplotypes < selected_haplotypes.len() {
        //     selected_haplotypes = selected_haplotypes[0..max_haplotypes].to_vec();
        // }
        // dbg!(&selected_haplotypes);

        dbg!(&extended_haplotypes);
        Ok(extended_haplotypes)
    }
    fn plot_prediction(
        &self,
        solution: &str,
        candidate_matrix_values: &Vec<(Vec<VariantStatus>, BitVec)>,
        haplotypes: &Vec<Haplotype>,
        variant_calls: &VariantCalls,
        best_variables: &Vec<f64>,
    ) -> Result<()> {
        let mut file_name = "".to_string();
        let json = include_str!("../../../templates/prediction.json");
        let mut blueprint: serde_json::Value = serde_json::from_str(json).unwrap();
        let mut plot_data_variants = Vec::new();
        let mut plot_data_haplotype_variants = Vec::new();
        let mut plot_data_haplotype_fractions = Vec::new();
        let mut plot_data_covered_variants = Vec::new();
        let mut plot_data_dataset_afd = Vec::new();

        if &solution == &"lp" {
            for ((genotype_matrix, coverage_matrix), (variant_id, (af, _))) in
                candidate_matrix_values.iter().zip(variant_calls.iter())
            {
                let mut counter = 0;
                for (i, _variable) in best_variables.iter().enumerate() {
                    if coverage_matrix[i as u64] {
                        counter += 1;
                    }
                }
                if counter == best_variables.len() {
                    for (i, (variable, haplotype)) in
                        best_variables.iter().zip(haplotypes.iter()).enumerate()
                    {
                        if genotype_matrix[i] == VariantStatus::Present {
                            plot_data_haplotype_fractions.push(DatasetHaplotypeFractions {
                                haplotype: haplotype.to_string(),
                                fraction: NotNan::new(*variable).unwrap(),
                            });
                            plot_data_haplotype_variants.push(DatasetHaplotypeVariants {
                                variant: *variant_id,
                                haplotype: haplotype.to_string(),
                            });
                            plot_data_variants.push(DatasetVariants {
                                variant: *variant_id,
                                vaf: af.clone(),
                            });
                        }
                    }
                }
            }
            file_name.push_str("lp_solution.json");
        } else if &solution == &"final" {
            candidate_matrix_values
                .iter()
                .zip(variant_calls.iter())
                .for_each(|((genotypes, covered), (variant_id, (af, afd)))| {
                    best_variables
                        .iter()
                        .zip(haplotypes.iter())
                        .enumerate()
                        .for_each(|(i, (fraction, haplotype))| {
                            if genotypes[i] == VariantStatus::Present && covered[i as u64] {
                                plot_data_haplotype_fractions.push(DatasetHaplotypeFractions {
                                    haplotype: haplotype.to_string(),
                                    fraction: NotNan::new(*fraction).unwrap(),
                                });
                                plot_data_haplotype_variants.push(DatasetHaplotypeVariants {
                                    variant: *variant_id,
                                    haplotype: haplotype.to_string(),
                                });
                                plot_data_variants.push(DatasetVariants {
                                    variant: *variant_id,
                                    vaf: *af,
                                });
                                //addition of one more rect plot for coverage matrix in addition to genotype matrix
                                //create the plot_data_covered_variants using only the variants that have GT:1 for at least one haplotype.
                                for (j, haplotype) in haplotypes.iter().enumerate() {
                                    if covered[j as u64] {
                                        plot_data_covered_variants.push(DatasetHaplotypeVariants {
                                            variant: *variant_id,
                                            haplotype: haplotype.to_string(),
                                        });
                                    }
                                }

                                //also add the heatmap for afd below the covered panels
                                for (allele_freq, prob) in afd.iter() {
                                    plot_data_dataset_afd.push(DatasetAfd {
                                        variant: *variant_id,
                                        allele_freq: *allele_freq,
                                        probability: prob.clone(),
                                    })
                                }
                            }
                        });
                });
            file_name.push_str("final_solution.json");
        }
        let plot_data_variants = json!(plot_data_variants);
        let plot_data_haplotype_variants = json!(plot_data_haplotype_variants);
        let plot_data_haplotype_fractions = json!(plot_data_haplotype_fractions);
        let plot_data_covered_variants = json!(plot_data_covered_variants);
        let plot_data_dataset_afd = json!(plot_data_dataset_afd);

        blueprint["datasets"]["variants"] = plot_data_variants;
        blueprint["datasets"]["haplotype_variants"] = plot_data_haplotype_variants;
        blueprint["datasets"]["haplotype_fractions"] = plot_data_haplotype_fractions;
        blueprint["datasets"]["covered_variants"] = plot_data_covered_variants;
        blueprint["datasets"]["allele_frequency_distribution"] = plot_data_dataset_afd;

        let mut parent = self.outcsv.clone().unwrap();
        parent.pop();
        fs::create_dir_all(&parent)?;
        let file = fs::File::create(parent.join(file_name)).unwrap();
        serde_json::to_writer(file, &blueprint)?;
        Ok(())
    }
    pub fn convert_to_g(&self) -> Result<BTreeMap<String, String>> {
        let mut reader = xml_reader::from_file(&self.xml)?;
        reader.trim_text(true);
        let mut buf = Vec::new();
        let mut alleles: Vec<String> = Vec::new();
        let mut confirmed: Vec<String> = Vec::new();
        let mut hla_g_groups: HashMap<i32, String> = HashMap::new(); //some hla alleles dont have g groups information in the xml file.
        let mut names_indices: Vec<i32> = Vec::new();
        let mut groups_indices: Vec<i32> = Vec::new();
        let mut counter = 0;
        loop {
            match reader.read_event_into(&mut buf) {
                Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
                Ok(Event::Eof) => break,
                Ok(Event::Start(e)) => match e.name().as_ref() {
                    b"allele" => {
                        alleles.push(
                            e.attributes()
                                .map(|a| String::from_utf8(a.unwrap().value.to_vec()))
                                .collect::<Vec<_>>()[1]
                                .as_ref()
                                .unwrap()
                                .split("-") //"HLA-" don't take the HLA prefix
                                .collect::<Vec<&str>>()[1]
                                .to_string(),
                        ); //allele_name is held in index 1, note: don't use expanded_name.
                        names_indices.push(counter.clone());
                        counter += 1;
                    }
                    _ => (),
                },
                Ok(Event::Empty(e)) => match e.name().as_ref() {
                    b"releaseversions" => confirmed.push(
                        e.attributes()
                            .map(|a| String::from_utf8(a.unwrap().value.to_vec()))
                            .collect::<Vec<_>>()[4]
                            .as_ref()
                            .unwrap()
                            .to_string(), //index 4 holds the Confirmed info
                    ),
                    b"hla_g_group" => {
                        groups_indices.push(counter.clone());
                        hla_g_groups.insert(
                            counter.clone(),
                            e.attributes()
                                .map(|a| String::from_utf8(a.unwrap().value.to_vec()))
                                .collect::<Vec<_>>()[0]
                                .as_ref()
                                .unwrap()
                                .to_string(), //index 0 holds the status info
                        );
                    }
                    _ => (),
                },
                _ => (),
            }
            // if we don't keep a borrow elsewhere, we can clear the buffer to keep memory usage low
            buf.clear();
        }
        assert_eq!(alleles.len(), confirmed.len());
        let mut filtered_alleles = Vec::new();
        let mut filtered_confirmed = Vec::new();
        hla_g_groups.iter().for_each(|(index, _)| {
            filtered_alleles.push(alleles[*index as usize - 1].clone());
            filtered_confirmed.push(confirmed[*index as usize - 1].clone());
        });
        assert_eq!(filtered_alleles.len(), filtered_confirmed.len());
        assert_eq!(filtered_alleles.len(), hla_g_groups.len());

        let mut g_to_alleles: BTreeMap<String, String> = BTreeMap::new();
        let g_names: Vec<String> = hla_g_groups.values().cloned().collect();
        let _unconfirmed_alleles = filtered_alleles
            .iter()
            .zip(filtered_confirmed.iter())
            .zip(g_names.iter())
            .filter(|((_allele, c), _g_group)| c == &"Confirmed")
            .for_each(|((allele, _c), g_group)| {
                g_to_alleles.insert(allele.clone(), g_group.to_string());
            });
        dbg!(&g_to_alleles);
        Ok(g_to_alleles)
    }
    fn plot_densities(
        &self,
        event_posteriors: &Vec<(HaplotypeFractions, LogProb)>,
        final_haplotypes: &Vec<Haplotype>,
        file_prefix: &str,
    ) -> Result<()> {
        let file_name = format!("{}_solutions.json", file_prefix.to_string());
        let json = include_str!("../../../templates/densities.json");
        let mut blueprint: serde_json::Value = serde_json::from_str(json).unwrap();
        let mut plot_data_fractions = Vec::new();
        let mut plot_density = Vec::new();

        //take first 10 solutions, if the events are less than that, then take the length
        let plot_first_events = 10;
        let num_events = event_posteriors.len();
        let mut new_event_posteriors = event_posteriors.clone();
        if num_events > plot_first_events {
            new_event_posteriors = new_event_posteriors[0..plot_first_events].to_vec();
        }
        for (i, (fractions, logprob)) in new_event_posteriors.iter().enumerate() {
            plot_density.push(DatasetDensitySolution {
                density: logprob.exp().clone(),
                solution_number: i,
            });
            for (f, h) in fractions.iter().zip(final_haplotypes.iter()) {
                //fractions and final haplotypes are already ordered.
                if f != &NotNan::new(0.0).unwrap() {
                    plot_data_fractions.push(DatasetHaplotypeFractionsWsolution {
                        haplotype: h.to_string(),
                        fraction: f.clone(),
                        solution_number: i,
                    });
                }
            }
        }

        //write to json
        let plot_data_fractions = json!(plot_data_fractions);
        let plot_density = json!(plot_density);
        blueprint["datasets"]["haplotype_fractions"] = plot_data_fractions;
        blueprint["datasets"]["densities"] = plot_density;

        let mut parent = self.outcsv.clone().unwrap();
        parent.pop();
        let file = fs::File::create(parent.join(file_name)).unwrap();
        serde_json::to_writer(file, &blueprint)?;
        Ok(())
    }
}

#[derive(Serialize, Debug)]
pub(crate) struct DatasetVariants {
    variant: VariantID,
    vaf: f32,
}
#[derive(Serialize, Debug)]
pub(crate) struct DatasetHaplotypeVariants {
    variant: VariantID,
    haplotype: String,
}
#[derive(Serialize, Debug)]
pub(crate) struct DatasetHaplotypeFractions {
    haplotype: String,
    fraction: AlleleFreq,
}

#[derive(Serialize, Debug)]
pub(crate) struct DatasetAfd {
    variant: VariantID,
    allele_freq: AlleleFreq,
    probability: f64,
}

#[derive(Serialize, Debug)]
pub(crate) struct DatasetHaplotypeFractionsWsolution {
    haplotype: String,
    fraction: AlleleFreq,
    solution_number: usize,
}

#[derive(Serialize, Debug)]
pub(crate) struct DatasetDensitySolution {
    density: f64,
    solution_number: usize,
}

// #[derive(Debug, Clone, Copy, Ord, PartialOrd, Eq, PartialEq)]
// pub(crate) struct KallistoEstimate {
//     pub count: NotNan<f64>,
//     pub dispersion: NotNan<f64>,
// }

// #[derive(Debug, Clone, Derefable, DerefMut)]
// pub(crate) struct KallistoEstimates(#[deref] BTreeMap<Haplotype, KallistoEstimate>);

// impl KallistoEstimates {
//     /// Generate new instance.
//     pub(crate) fn new(
//         hdf5_reader: &hdf5::File,
//         min_norm_counts: f64,
//         haplotypes: &Vec<Haplotype>,
//     ) -> Result<Self> {
//         let seqnames = Self::filter_seqnames(hdf5_reader, min_norm_counts)?;
//         let ids = hdf5_reader
//             .dataset("aux/ids")?
//             .read_1d::<hdf5::types::FixedAscii<255>>()?;
//         let num_bootstraps = hdf5_reader.dataset("aux/num_bootstrap")?.read_1d::<i32>()?;
//         let seq_length = hdf5_reader.dataset("aux/lengths")?.read_1d::<f64>()?;
//         let mut estimates = BTreeMap::new();
//         for seqname in seqnames {
//             if haplotypes.contains(&Haplotype(seqname.clone())) {
//                 let index = ids.iter().position(|x| x.as_str() == seqname).unwrap();
//                 let mut bootstraps = Vec::new();
//                 for i in 0..num_bootstraps[0] {
//                     let dataset = hdf5_reader.dataset(&format!("bootstrap/bs{i}", i = i))?;
//                     let est_counts = dataset.read_1d::<f64>()?;
//                     let norm_counts = est_counts / &seq_length;
//                     let norm_counts = norm_counts[index];
//                     bootstraps.push(norm_counts);
//                 }

//                 //mean
//                 let sum = bootstraps.iter().sum::<f64>();
//                 let count = bootstraps.len();
//                 let m = sum / count as f64;

//                 //std dev
//                 let variance = bootstraps
//                     .iter()
//                     .map(|value| {
//                         let diff = m - (*value as f64);
//                         diff * diff
//                     })
//                     .sum::<f64>()
//                     / count as f64;
//                 let std = variance.sqrt();
//                 let t = std / m;
//                 //retrieval of mle
//                 let mle_dataset = hdf5_reader.dataset("est_counts")?.read_1d::<f64>()?;
//                 let mle_norm = mle_dataset / &seq_length; //normalized mle counts by length
//                 let m = mle_norm[index];
//                 estimates.insert(
//                     Haplotype(seqname.clone()),
//                     KallistoEstimate {
//                         dispersion: NotNan::new(t).unwrap(),
//                         count: NotNan::new(m).unwrap(),
//                     },
//                 );
//             }
//         }
//         Ok(KallistoEstimates(estimates))
//     }

//     //Return top N estimates according to --max-haplotypes
//     fn select_haplotypes(self, max_haplotypes: i64) -> Result<Self> {
//         let mut kallisto_estimates = self.clone();
//         // kallisto_estimates.retain(|k, _| haplotypes.contains(&k));
//         let mut estimates_vec: Vec<(&Haplotype, &KallistoEstimate)> =
//             kallisto_estimates.iter().collect();
//         estimates_vec.sort_by(|a, b| b.1.count.cmp(&a.1.count));
//         if estimates_vec.len() >= max_haplotypes.try_into().unwrap() {
//             let topn = estimates_vec[0..max_haplotypes as usize].to_vec();
//             let mut top_estimates = BTreeMap::new();
//             for (key, value) in topn {
//                 top_estimates.insert(key.clone(), *value);
//             }
//             Ok(KallistoEstimates(top_estimates))
//         } else {
//             Ok(self)
//         }
//     }
//     //Return a vector of filtered seqnames according to --min-norm-counts.
//     fn filter_seqnames(hdf5_reader: &hdf5::File, min_norm_counts: f64) -> Result<Vec<String>> {
//         let ids = hdf5_reader
//             .dataset("aux/ids")?
//             .read_1d::<hdf5::types::FixedAscii<255>>()?;
//         let est_counts = hdf5_reader.dataset("est_counts")?.read_1d::<f64>()?;
//         let seq_length = hdf5_reader.dataset("aux/lengths")?.read_1d::<f64>()?; //these two variables arrays have the same length.
//         let norm_counts = est_counts / seq_length;
//         let mut filtered_haplotypes: Vec<String> = Vec::new();
//         for (num, id) in norm_counts.iter().zip(ids.iter()) {
//             if num > &min_norm_counts {
//                 filtered_haplotypes.push(id.to_string());
//             }
//         }
//         Ok(filtered_haplotypes)
//     }
// }

fn collect_constraints_and_variants(
    candidate_matrix: &CandidateMatrix,
    haplotypes: &Vec<Haplotype>,
    variant_calls: &VariantCalls,
    variables: &Vec<Variable>,
    constraints: &mut Vec<Expression>,
) -> Result<HashMap<Haplotype, Vec<VariantID>>> {
    let candidate_matrix_values: Vec<(Vec<VariantStatus>, BitVec)> =
        candidate_matrix.values().cloned().collect();
    //collect haplotype-to-variants information
    let mut haplotype_dict: HashMap<Haplotype, Vec<VariantID>> =
        haplotypes.iter().map(|h| (h.clone(), vec![])).collect();
    //variant-wise iteration
    let mut expr = Expression::from_other_affine(0.); // A constant expression
    for ((genotype_matrix, coverage_matrix), (variant, (af, _))) in
        candidate_matrix_values.iter().zip(variant_calls.iter())
    {
        let mut fraction_cont = Expression::from_other_affine(0.);
        let _prime_fraction_cont = Expression::from_other_affine(0.);
        let _vaf = Expression::from_other_affine(0.);
        let mut counter = 0;
        for (i, _variable) in variables.iter().enumerate() {
            if coverage_matrix[i as u64] {
                counter += 1;
            }
        }
        if counter == variables.len() {
            for (i, (variable, haplotype)) in variables.iter().zip(haplotypes.iter()).enumerate() {
                if genotype_matrix[i] == VariantStatus::Present {
                    fraction_cont += *variable;
                    let mut existing = haplotype_dict.get(&haplotype).unwrap().clone();
                    existing.push(variant.clone());
                    haplotype_dict.insert(haplotype.clone(), existing);
                }
            }
            let expr_to_add = fraction_cont - af.clone().into_expression();
            constraints.push(expr_to_add.clone());
            expr += expr_to_add;
        }
    }
    Ok(haplotype_dict)
}

//convert_to_two_field function converts the event posteriors that contain three-field info by default, to two-field information
//by summing densities of events that have identical explanation with the first two fields
fn convert_to_two_field(
    event_posteriors: &Vec<(HaplotypeFractions, LogProb)>,
    haplotypes: &Vec<Haplotype>,
) -> Result<(Vec<Haplotype>, Vec<(HaplotypeFractions, LogProb)>)> {
    // dbg!(&event_posteriors);
    let mut event_posteriors_map: Vec<(BTreeMap<Haplotype, NotNan<f64>>, LogProb)> = Vec::new();
    dbg!(&event_posteriors);
    for (fractions, logprob) in event_posteriors.iter() {
        //firstly, initiate a map for haplotype and fraction info for each event
        //by having zero fraction as first values
        let mut haplotype_to_fraction_new: BTreeMap<Haplotype, NotNan<f64>> = haplotypes
            .iter()
            .map(|h| {
                let splitted: Vec<&str> = h.split(':').collect();
                let two_field = format!("{}:{}", splitted[0].to_string(), splitted[1]);
                (Haplotype(two_field.clone()), NotNan::new(0.00).unwrap())
            })
            .collect();

        // secondly, enter the fraction per haplotype to fraction map
        // if the value has not been updated before to a value that is greater than 0.0
        // that is done to prevent overwriting the haplotype that has a fraction greater
        // than 0 to 0 because of having the identical first two fields
        for (fraction, haplotype) in fractions.iter().zip(haplotypes.iter()) {
            let splitted: Vec<&str> = haplotype.split(':').collect();
            let two_field = format!("{}:{}", splitted[0].to_string(), splitted[1]);
            let two_field = Haplotype(two_field);
            if haplotype_to_fraction_new[&two_field] == NotNan::new(0.00).unwrap() {
                haplotype_to_fraction_new.insert(two_field.clone(), *fraction);
            } else {
                // this is to ensure that haplotypes with identical two_fields do not have separate records
                let sum_of_two = haplotype_to_fraction_new[&two_field].clone() + fraction.clone();
                haplotype_to_fraction_new.insert(two_field.clone(), sum_of_two);
            }
        }
        event_posteriors_map.push((haplotype_to_fraction_new.clone(), *logprob));
    }

    //last, create a map for haplotype-fractions to logprob,
    //in order to sum all logprobs belonging to same haplotype-fractions
    let mut hf_to_logprob: BTreeMap<BTreeMap<Haplotype, NotNan<f64>>, LogProb> = BTreeMap::new();
    // dbg!(&event_posteriors_map);
    for (hf, logprob) in event_posteriors_map.iter() {
        if hf_to_logprob.contains_key(&hf) {
            let new_logprob = LogProb::ln_sum_exp(&vec![hf_to_logprob[&hf].clone(), *logprob]);
            hf_to_logprob.insert(hf.clone(), new_logprob.clone());
        } else {
            hf_to_logprob.insert(hf.clone(), logprob.clone());
        }
    }
    // dbg!(&hf_to_logprob);
    //logprob doesn't implement Ord, so, convert the map to a vector of tuples starting with logprob
    let mut logprob_and_hf = Vec::new();
    for (hf, logprob) in hf_to_logprob.iter() {
        logprob_and_hf.push((logprob, hf));
    }
    logprob_and_hf.sort_by(|a, b| b.partial_cmp(a).unwrap_or(Ordering::Equal));

    //convert the final construct to the same type with the input of the function
    //(event_posteriors) and finally return new haplotypes with two-field information
    //in addition to event posteriors
    let (_lp, map) = &logprob_and_hf[0];
    let final_haplotypes: Vec<Haplotype> = map.keys().cloned().collect();
    let mut event_posteriors_two_field = Vec::new();
    for (lp, map) in logprob_and_hf.iter() {
        let haplotype_fractions = HaplotypeFractions(map.values().cloned().collect());
        event_posteriors_two_field.push((haplotype_fractions, **lp));
    }
    dbg!(&event_posteriors_two_field);
    Ok((final_haplotypes, event_posteriors_two_field))
    // Ok(())
}