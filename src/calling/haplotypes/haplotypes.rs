use crate::model::{AlleleFreq, Data, HaplotypeFractions};
use anyhow::Result;
use bio::stats::{probs::LogProb, PHREDProb, Prob};
use bv::BitVec;

use derefable::Derefable;

use derive_deref::DerefMut;

use ordered_float::NotNan;
use ordered_float::OrderedFloat;

use rust_htslib::bcf::{
    self,
    record::GenotypeAllele::{Phased, Unphased},
    Read,
};

use good_lp::IntoAffineExpression;
use good_lp::*;
use good_lp::{variable, Expression};

use serde::Serialize;
use serde_json::json;
use std::collections::{BTreeMap, HashMap};

use std::fs;
use std::str::FromStr;
use std::{path::PathBuf, str};

#[derive(Derefable, Debug, Copy, Clone, PartialEq, Eq, Hash, Ord, PartialOrd, Serialize)]
pub struct VariantID(#[deref] pub i32);

#[derive(Derefable, Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize)]
pub struct Haplotype(#[deref] pub String);

#[derive(Debug, Clone, Derefable)]
pub struct AlleleFreqDist(#[deref] BTreeMap<AlleleFreq, f64>);

impl AlleleFreqDist {
    pub fn vaf_query(&self, vaf: &AlleleFreq) -> Option<LogProb> {
        if self.contains_key(&vaf) {
            Some(LogProb::from(PHREDProb(*self.get(&vaf).unwrap())))
        } else {
            let (x_0, y_0) = self.range(..vaf).next_back().unwrap();
            let (x_1, y_1) = self.range(vaf..).next().unwrap();
            let density =
                NotNan::new(*y_0).unwrap() + (*vaf - *x_0) * (*y_1 - *y_0) / (*x_1 - *x_0); //calculation of density for given vaf by linear interpolation
            Some(LogProb::from(PHREDProb(NotNan::into_inner(density))))
        }
    }
}

#[derive(Derefable, Debug, Clone)]
pub struct CandidateMatrix(#[deref] BTreeMap<VariantID, (Vec<VariantStatus>, BitVec)>);

impl CandidateMatrix {
    pub fn new(
        haplotype_variants: &HaplotypeVariants,
        // haplotypes: &Vec<Haplotype>,
    ) -> Result<Self> {
        let mut candidate_matrix = BTreeMap::new();
        haplotype_variants.iter().for_each(|(variant_id, bmap)| {
            let mut haplotype_variants_gt = Vec::new();
            let mut haplotype_variants_c = BitVec::new();
            bmap.iter().for_each(|(_haplotype, (gt, c))| {
                haplotype_variants_c.push(*c);
                haplotype_variants_gt.push(gt.clone());
            });
            candidate_matrix.insert(*variant_id, (haplotype_variants_gt, haplotype_variants_c));
        });
        Ok(CandidateMatrix(candidate_matrix))
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum PriorTypes {
    Diploid,
    DiploidSubclonal,
    Uniform,
}

impl FromStr for PriorTypes {
    type Err = ();

    fn from_str(input: &str) -> Result<PriorTypes, Self::Err> {
        match input {
            "uniform" => Ok(PriorTypes::Uniform),
            "diploid" => Ok(PriorTypes::Diploid),
            "diploid-subclonal" => Ok(PriorTypes::DiploidSubclonal),
            _ => Err(()),
        }
    }
}

#[derive(Derefable, DerefMut, Debug, Clone)]
pub struct VariantCalls(#[deref] BTreeMap<VariantID, (f32, AlleleFreqDist)>); //The place of f32 is maximum a posteriori estimate of AF.

impl VariantCalls {
    pub fn new(variant_calls: &mut bcf::Reader) -> Result<Self> {
        let mut calls = BTreeMap::new();
        for record_result in variant_calls.records() {
            let mut record = record_result?;
            record.unpack();
            let prob_absent = record.info(b"PROB_ABSENT").float().unwrap().unwrap()[0];
            let prob_absent_prob = Prob::from(PHREDProb(prob_absent.into()));
            let afd_utf = record.format(b"AFD").string()?;
            let afd = std::str::from_utf8(afd_utf[0]).unwrap();
            let read_depths = record.format(b"DP").integer().unwrap();
            if read_depths[0] != &[0]
            // && (&prob_absent_prob <= &Prob(0.05) || &prob_absent_prob >= &Prob(0.95))
            {
                //because some afd strings are just "." and that throws an error while splitting below.
                let variant_id: i32 = String::from_utf8(record.id())?.parse().unwrap();
                let af = (&*record.format(b"AF").float().unwrap()[0]).to_vec()[0];
                //dbg!(&af);
                let mut vaf_density = BTreeMap::new();
                for pair in afd.split(',') {
                    if let Some((vaf, density)) = pair.split_once("=") {
                        let (vaf, density): (AlleleFreq, f64) =
                            (vaf.parse().unwrap(), density.parse().unwrap());
                        vaf_density.insert(vaf, density);
                    }
                }
                calls.insert(VariantID(variant_id), (af, AlleleFreqDist(vaf_density)));
            }
        }
        Ok(VariantCalls(calls))
    }
    pub fn filter_variant_calls(&self, variants: &Vec<VariantID>) -> Result<Self> {
        let mut variant_calls_filtered = self.clone();
        for (v, _) in self.iter() {
            if !variants.contains(&v) {
                variant_calls_filtered.remove_entry(&v);
            }
        }
        Ok(variant_calls_filtered)
    }
}

#[derive(Clone, Debug, Eq, PartialEq, PartialOrd)]
pub enum VariantStatus {
    Present,
    NotPresent,
    Unknown,
}

#[derive(Derefable, Debug, Clone, PartialEq, Eq, PartialOrd, DerefMut)]
pub struct HaplotypeVariants(
    #[deref] BTreeMap<VariantID, BTreeMap<Haplotype, (VariantStatus, bool)>>,
);

impl HaplotypeVariants {
    pub fn new(
        //observations: &mut bcf::Reader,
        haplotype_variants: &mut bcf::Reader,
        filtered_ids: &Vec<VariantID>,
        //max_haplotypes: &usize,
    ) -> Result<Self> {
        let mut variant_records = BTreeMap::new();
        for record_result in haplotype_variants.records() {
            let record = record_result?;
            let variant_id: VariantID = VariantID(String::from_utf8(record.id())?.parse().unwrap());
            if filtered_ids.contains(&variant_id) {
                let header = record.header();
                let gts = record.genotypes()?;
                let loci = record.format(b"C").integer().unwrap();
                let mut matrices = BTreeMap::new();
                for (index, haplotype) in header.samples().iter().enumerate() {
                    let haplotype = Haplotype(str::from_utf8(haplotype).unwrap().to_string());
                    //generate phased genotypes.
                    for gta in gts.get(index).iter().skip(1) {
                        //maternal and paternal gts will be the same in the vcf i.e. 0|0 and 1|1
                        if *gta == Unphased(1) || *gta == Phased(1) {
                            matrices.insert(
                                haplotype.clone(),
                                (VariantStatus::Present, loci[index] == &[1]),
                            );
                        } else {
                            matrices.insert(
                                haplotype.clone(),
                                (VariantStatus::NotPresent, loci[index] == &[1]),
                            );
                        }
                    }
                }
                variant_records.insert(variant_id, matrices);
            }
        }
        Ok(HaplotypeVariants(variant_records))
    }

    pub fn find_plausible_haplotypes(
        &self,
        _variant_calls: &VariantCalls,
        haplotypes: &Vec<Haplotype>,
    ) -> Result<Self> {
        let mut new_haplotype_variants: BTreeMap<
            VariantID,
            BTreeMap<Haplotype, (VariantStatus, bool)>,
        > = BTreeMap::new();
        for (variant, matrix_map) in self.iter() {
            let mut new_matrix_map = BTreeMap::new();
            for (haplotype_m, (variant_status, coverage_status)) in matrix_map {
                if haplotypes.contains(&haplotype_m) {
                    //fix: filter for haplotypes in the haplotypes list
                    new_matrix_map.insert(
                        haplotype_m.clone(),
                        (variant_status.clone(), coverage_status.clone()),
                    );
                }
            }
            new_haplotype_variants.insert(variant.clone(), new_matrix_map);
        }
        Ok(HaplotypeVariants(new_haplotype_variants))
    }

    pub fn find_common_variants(
        &self,
        variant_calls: &VariantCalls,
        haplotypes: &Vec<Haplotype>,
    ) -> Result<Vec<VariantID>> {
        let candidate_matrix_values: Vec<(Vec<VariantStatus>, BitVec)> = CandidateMatrix::new(self)
            .unwrap()
            .values()
            .cloned()
            .collect();
        let mut common_variants = Vec::new();
        for ((_genotype_matrix, coverage_matrix), (variant, (_af, _))) in
            candidate_matrix_values.iter().zip(variant_calls.iter())
        {
            let mut counter = 0;
            for (i, _haplotype) in haplotypes.iter().enumerate() {
                if coverage_matrix[i as u64] {
                    counter += 1;
                }
            }
            if counter == haplotypes.len() {
                common_variants.push(variant.clone());
            }
        }
        Ok(common_variants)
    }
    pub fn filter_haplotype_variants(&self, variants: &Vec<VariantID>) -> Result<Self> {
        let mut haplotype_variants_filtered = self.clone();
        for (v, _) in self.iter() {
            if !variants.contains(&v) {
                haplotype_variants_filtered.remove_entry(&v);
            }
        }
        Ok(haplotype_variants_filtered)
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

pub fn plot_prediction(
    outdir: &PathBuf,
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

    let mut parent = outdir.clone();
    parent.pop();
    fs::create_dir_all(&parent)?;
    let file = fs::File::create(parent.join(file_name)).unwrap();
    serde_json::to_writer(file, &blueprint)?;
    Ok(())
}

pub fn linear_program(
    outdir: &PathBuf,
    candidate_matrix: &CandidateMatrix,
    haplotypes: &Vec<Haplotype>,
    variant_calls: &VariantCalls,
    lp_cutoff: f64,
    num_variant_distance: i64,
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
    let t_vars: Vec<Variable> = problem.add_vector(variable().min(0.0).max(1.0), constraints.len());

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
    plot_prediction(
        outdir,
        &"lp",
        &candidate_matrix_values,
        &haplotypes,
        &variant_calls,
        &best_variables,
    )?;

    //extend haplotypes found by linear program, add haplotypes that have the same variants to the final list.
    //then sort by hamming distance, take the closest x additional alleles according to 'num_variant_distance'.
    //this is done by storing only the variants that have GT:1 and C:1 for all haplotypes in haplotype_dict and remaining variants are not included.

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
                    let mut difference = vec![];
                    for i in haplotype_variants.iter() {
                        if !variants.contains(&i) {
                            difference.push(i);
                        }
                    }
                    if (difference.len() as i64 <= num_variant_distance)
                        && ((variants.len() as i64 - haplotype_variants.len() as i64).abs()
                            <= num_variant_distance)
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

pub fn write_results(
    outdir: &PathBuf,
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
                    vaf_sum = NotNan::new((vaf_sum * NotNan::new(100.0).unwrap()).round()).unwrap()
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
    let mut wtr = csv::Writer::from_path(outdir)?;
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

pub fn collect_constraints_and_variants(
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

pub fn plot_densities(
    outdir: &PathBuf,
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

    let mut parent = outdir.clone();
    parent.pop();
    let file = fs::File::create(parent.join(file_name)).unwrap();
    serde_json::to_writer(file, &blueprint)?;
    Ok(())
}
