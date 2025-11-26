use crate::calling::haplotypes::haplotypes;
use crate::calling::haplotypes::haplotypes::filter_variants_for_best_solution_plot;
use crate::calling::haplotypes::haplotypes::get_event_posteriors;
use crate::calling::haplotypes::haplotypes::output_empty_output;
use crate::calling::haplotypes::haplotypes::{
    CandidateMatrix, Haplotype, HaplotypeVariants, VariantCalls,
};
use crate::model::HaplotypeFractions;
use anyhow::Result;
use bio::stats::probs::LogProb;

use core::cmp::Ordering;

use derive_builder::Builder;

use ordered_float::NotNan;

use quick_xml::events::Event;
use quick_xml::reader::Reader as xml_reader;

use rust_htslib::bcf::{self};

use std::collections::{BTreeMap, HashMap};

use std::fs;

use std::{path::PathBuf, str};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    haplotype_variants: bcf::Reader,
    variant_calls: bcf::Reader,
    xml: PathBuf,
    output_folder: PathBuf,
    prior: String,
    // common_variants: bool,
    lp_cutoff: f64,
    enable_equivalence_class_constraint: bool,
    extend_haplotypes: Option<bool>,
    threshold_equivalence_class: usize,
    num_extend_haplotypes: i64,
    num_constraint_haplotypes: i32,
    output_lp_datavzrd: bool,
    sample_name: Option<String>,
    enforce_given_alleles: Option<Vec<String>>,
}

impl Caller {
    pub fn call(&mut self) -> Result<()> {
        //Step 1: Prepare data and compute the model
        //initially prepare haplotype_variants and variant_calls
        let variant_calls = VariantCalls::new(&mut self.variant_calls, &self.sample_name)?;

        //write blank plots and tsv table if no variants are available.
        if variant_calls.len() == 0 {
            output_empty_output(&self.output_folder).unwrap();
            Ok(())
        } else {
            let mut haplotype_variants = HaplotypeVariants::new(&mut self.haplotype_variants)?;

            //filter candidates vcf based on optional given input set of alleles (3-field-resolution)
            if let Some(input_alleles) = &self.enforce_given_alleles {
                // dbg!(&input_alleles);
                haplotype_variants =
                    haplotype_variants.filter_for_haplotype_prefixes(&input_alleles)?;
                // dbg!(&haplotype_variants);
            }

            let (event_posteriors, all_haplotypes, data) = get_event_posteriors(
                &self.output_lp_datavzrd,
                &haplotype_variants,
                variant_calls,
                &"hla",
                &self.prior,
                &self.output_folder,
                self.extend_haplotypes.unwrap_or(true),
                self.num_extend_haplotypes,
                self.num_constraint_haplotypes,
                self.lp_cutoff,
                self.enable_equivalence_class_constraint,
                Some(self.threshold_equivalence_class),
            )?;
            // dbg!(&event_posteriors, &all_haplotypes);

            //plot the best solution as final solution plot
            let best_fractions = event_posteriors
                .iter()
                .next()
                .unwrap()
                .0
                .iter()
                .map(|f| NotNan::into_inner(*f))
                .collect::<Vec<f64>>();

            //output best solution plot for only nonzero fraction haplotypes
            let nonzero_haplotype_fractions: BTreeMap<Haplotype, f64> =
                get_nonzero_haplotype_fractions(&all_haplotypes, &best_fractions);

            //collect haplotype names and fractions separately to be used later twice
            let filtered_haplotypes = nonzero_haplotype_fractions.keys().cloned().collect();
            let filtered_fractions = nonzero_haplotype_fractions.values().cloned().collect();
            // dbg!(&filtered_haplotypes, &filtered_fractions);

            //filter candidate matrix based on nonzero haplotype fractions
            let filtered_candidate_matrix = CandidateMatrix::new(
                &haplotype_variants
                    .filter_for_haplotypes(&filtered_haplotypes)
                    .unwrap(),
            )
            .unwrap();

            //plot best solution plot and only display variants that are found in one of the haplotypes.
            let (best_solution_matrix, best_solution_variant_calls) =
                filter_variants_for_best_solution_plot(
                    &filtered_candidate_matrix,
                    &data.variant_calls,
                );

            haplotypes::plot_prediction(
                &self.output_lp_datavzrd,
                &self.output_folder,
                &"final",
                &best_solution_matrix,
                &filtered_haplotypes,
                &best_solution_variant_calls,
                &filtered_fractions,
            )?;

            //write results to tsv
            haplotypes::write_results(
                &self.output_folder.join(&"predictions.csv"),
                &data,
                &event_posteriors,
                &all_haplotypes,
                self.prior.clone(),
                false,
            )?;

            // arrow plot

            //filter out the variants that are not within the range of the locus.
            let first_haplotype = nonzero_haplotype_fractions.keys().next().unwrap();
            let (locus_start, locus_end) = first_haplotype.get_coordinates_for_haplotype();
            let (variant_calls_in_locus, candidate_matrix_in_locus) = &data
                .variant_calls
                .filter_variants_in_range(&filtered_candidate_matrix, locus_start, locus_end)?;

            //the haplotype order is preserved in the keys of nonzero_haplotype_fractions
            haplotypes::get_arrow_plot(
                &self.output_folder,
                candidate_matrix_in_locus,
                &nonzero_haplotype_fractions,
                variant_calls_in_locus,
            );

            //second: 2-field
            let (two_field_haplotypes, two_field_event_posteriors) =
                convert_to_two_field(&event_posteriors, &all_haplotypes)?;
            let mut path_for_two_fields = PathBuf::from(&self.output_folder);
            path_for_two_fields.push("2-field.csv");
            haplotypes::write_results(
                &path_for_two_fields,
                &data,
                &two_field_event_posteriors,
                &two_field_haplotypes,
                self.prior.clone(),
                false,
            )?;

            //plot first 10 posteriors of orthanq output
            haplotypes::plot_densities(
                &self.output_folder,
                &event_posteriors,
                &all_haplotypes,
                "3_field",
            )?;
            haplotypes::plot_densities(
                &self.output_folder,
                &two_field_event_posteriors,
                &two_field_haplotypes,
                "2_field",
            )?;

            //write table for G groups of HLA alleles, for HLA alleles with None G group in the XML table, we write the haplotype name back.
            //as a hint, successfuly converted G groups will have G in the end, while the ones with no G group will not have one.
            let mut converted_name = PathBuf::from(&self.output_folder);
            converted_name.push("G_groups.csv");
            let allele_to_g_groups = self.convert_to_g().unwrap();
            // dbg!(&allele_to_g_groups);
            // dbg!(&all_haplotypes);
            let mut final_haplotypes_converted: Vec<Haplotype> = Vec::new();

            //the haplotype can either be found with the same name in the xml file or it can start with it.
            //this is because we only use 3-field resolution of the haplotypes. This means, some haplotypes can
            // directly match, e.g. 24:436 while some e.g. A*24:03:01 have longer names in the xml.
            for haplotype in &all_haplotypes {
                let mut found_match = false;

                for (allele, g_group) in &allele_to_g_groups {
                    if allele == &haplotype.to_string()
                        || allele.starts_with(&haplotype.to_string())
                    {
                        if g_group == "None" {
                            //in case the allele has "None" in the g group field in the xml.
                            final_haplotypes_converted.push(haplotype.clone());
                        } else {
                            final_haplotypes_converted.push(Haplotype(g_group.to_string()));
                        }
                        found_match = true;
                        break;
                    }
                }

                if !found_match {
                    //in case the allele does not have a corresponding g group field in the xml.
                    final_haplotypes_converted.push(Haplotype(haplotype.to_string()));
                }
            }

            // dbg!(&final_haplotypes_converted);
            haplotypes::write_results(
                &converted_name,
                &data,
                &event_posteriors,
                &final_haplotypes_converted,
                self.prior.clone(),
                false,
            )?;
            Ok(())
        }
    }

    pub fn convert_to_g(&self) -> Result<BTreeMap<String, String>> {
        let mut reader = xml_reader::from_file(&self.xml)?;
        reader.trim_text(true);
        let mut buf = Vec::new();
        let mut allele_names: Vec<String> = Vec::new();
        let mut confirmed: Vec<String> = Vec::new();
        let mut hla_g_groups: HashMap<i32, String> = HashMap::new(); //some hla alleles dont have g groups information in the xml file.
        let mut groups_indices: Vec<i32> = Vec::new();
        //we keep track of each allele-g group pair with a counter to be used as index
        let mut counter = 0;
        loop {
            match reader.read_event_into(&mut buf) {
                Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
                Ok(Event::Eof) => break,
                Ok(Event::Start(e)) => match e.name().as_ref() {
                    b"allele" => {
                        let mut id_value: Option<String> = None;
                        let mut name_value: Option<String> = None;

                        for attr in e.attributes().flatten() {
                            if let Ok(key) = std::str::from_utf8(attr.key.as_ref()) {
                                if let Ok(val) = std::str::from_utf8(&attr.value) {
                                    match key {
                                        "id" => id_value = Some(val.to_string()),
                                        "name" => name_value = Some(val.to_string()),
                                        _ => {}
                                    }
                                }
                            }
                        }
                        match (id_value, name_value) {
                            (Some(id), Some(name)) => {
                                //clean up the allele name by removing the "HLA-" prefix if present
                                let cleaned_name = if name.contains('-') {
                                    name.split('-').nth(1).unwrap_or(&name).to_string()
                                } else {
                                    name
                                };
                                allele_names.push(cleaned_name);

                                counter += 1;
                            }
                            (id_opt, name_opt) => {
                                eprintln!(
                                    "Warning: missing attribute{}{} in <allele> element",
                                    if id_opt.is_none() { " 'id'" } else { "" },
                                    if name_opt.is_none() { " 'name'" } else { "" }
                                );
                            }
                        }
                    }
                    _ => (),
                },
                Ok(Event::Empty(e)) => match e.name().as_ref() {
                    b"hla_g_group" => {
                        let mut status_value: Option<String> = None;

                        for attr in e.attributes().flatten() {
                            if let Ok(key) = std::str::from_utf8(attr.key.as_ref()) {
                                if key == "status" {
                                    if let Ok(val) = std::str::from_utf8(&attr.value) {
                                        status_value = Some(val.to_string());
                                    }
                                }
                            }
                        }

                        if let Some(status) = status_value {
                            hla_g_groups.insert(counter, status);
                        } else {
                            eprintln!(
                                "Warning: No 'status' attribute found for hla_g_group at index {}",
                                counter
                            );
                        }
                    }
                    _ => (),
                },
                _ => (),
            }
            buf.clear();
        }

        let mut allele_to_g: BTreeMap<String, String> = BTreeMap::new();
        for (idx, g_group) in &hla_g_groups {
            if let Some(allele) = allele_names.get((*idx as usize) - 1) {
                allele_to_g.insert(allele.clone(), g_group.clone());
            }
        }
        // dbg!(&allele_to_g)
        Ok(allele_to_g)
    }
}

//convert_to_two_field function converts the event posteriors that contain three-field info by default, to two-field information
//by summing densities of events that have identical explanation with the first two fields
fn convert_to_two_field(
    event_posteriors: &Vec<(HaplotypeFractions, LogProb)>,
    haplotypes: &Vec<Haplotype>,
) -> Result<(Vec<Haplotype>, Vec<(HaplotypeFractions, LogProb)>)> {
    let mut event_posteriors_map: Vec<(BTreeMap<Haplotype, NotNan<f64>>, LogProb)> = Vec::new();
    // dbg!(&event_posteriors);
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
    for (hf, logprob) in event_posteriors_map.iter() {
        if hf_to_logprob.contains_key(&hf) {
            let new_logprob = LogProb::ln_sum_exp(&vec![hf_to_logprob[&hf].clone(), *logprob]);
            hf_to_logprob.insert(hf.clone(), new_logprob.clone());
        } else {
            hf_to_logprob.insert(hf.clone(), logprob.clone());
        }
    }
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
    // dbg!(&event_posteriors_two_field);
    Ok((final_haplotypes, event_posteriors_two_field))
    // Ok(())
}

fn get_nonzero_haplotype_fractions(
    haplotypes: &[Haplotype],
    fractions: &[f64],
) -> BTreeMap<Haplotype, f64> {
    haplotypes
        .iter()
        .zip(fractions.iter())
        .filter_map(|(hap, freq)| {
            if *freq != 0.0 {
                Some((hap.clone(), *freq))
            } else {
                None
            }
        })
        .collect()
}
