use crate::calling::haplotypes::haplotypes;
use crate::calling::haplotypes::haplotypes::get_event_posteriors;
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
    outcsv: PathBuf,
    prior: String,
    // common_variants: bool,
    lp_cutoff: f64,
    enable_equivalence_class_constraint: bool,
    extend_haplotypes: bool,
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
            let mut parent = self.outcsv.clone();
            parent.pop();
            fs::create_dir_all(&parent)?;

            //write blank plots, required for the workflow!
            for file_name in vec![
                "lp_solution.json".to_string(),
                "final_solution.json".to_string(),
                "2_field_solutions.json".to_string(),
                "3_field_solutions.json".to_string(),
            ] {
                let json = include_str!("../../../templates/final_prediction.json");
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
            let haplotype_variants = HaplotypeVariants::new(&mut self.haplotype_variants)?;
            let (event_posteriors, all_haplotypes, data) = get_event_posteriors(
                &haplotype_variants,
                variant_calls,
                &"hla",
                &self.prior,
                &self.outcsv,
                self.extend_haplotypes,
                self.num_extend_haplotypes,
                self.num_constraint_haplotypes,
                self.lp_cutoff,
                self.enable_equivalence_class_constraint,
                Some(self.threshold_equivalence_class),
            )?;
            // dbg!(&event_posteriors, &all_haplotypes);

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

            //second: 2-field
            let (two_field_haplotypes, two_field_event_posteriors) =
                convert_to_two_field(&event_posteriors, &all_haplotypes)?;
            let mut path_for_two_fields = PathBuf::from(&self.outcsv.parent().unwrap());
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
                &self.outcsv,
                &event_posteriors,
                &all_haplotypes,
                "3_field",
            )?;
            haplotypes::plot_densities(
                &self.outcsv,
                &two_field_event_posteriors,
                &two_field_haplotypes,
                "2_field",
            )?;

            //second: convert to G groups and write another table
            let mut converted_name = PathBuf::from(&self.outcsv.parent().unwrap());
            converted_name.push("G_groups.csv");
            let allele_to_g_groups = self.convert_to_g().unwrap();
            let mut final_haplotypes_converted: Vec<Haplotype> = Vec::new();
            all_haplotypes.iter().for_each(|haplotype| {
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
            //todo: retest this part and convert to true
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
                        let allele_name = e
                            .attributes()
                            .map(|a| String::from_utf8(a.unwrap().value.to_vec()))
                            .collect::<Vec<_>>()[1]
                            .clone()
                            .unwrap();
                        let mut allele_name_push = "".to_string();
                        //allele names starting with HLA contain '-' however, some aleles e.g. MICA do not contain it.
                        if allele_name.contains(&"-") {
                            allele_name_push =
                                allele_name.split("-").collect::<Vec<&str>>()[1].to_string()
                        } else {
                            allele_name_push = allele_name.clone()
                        }
                        alleles.push(allele_name_push); //allele_name is held in index 1, note: don't use expanded_name.
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
        Ok(g_to_alleles)
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
