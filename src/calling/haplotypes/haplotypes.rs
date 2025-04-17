use crate::model::{AlleleFreq, Data, HaplotypeFractions};
use crate::model::{Likelihood, Marginal, Posterior, Prior};
use anyhow::Result;
use bio::stats::bayesian::model::Model;
use bio::stats::{probs::LogProb, PHREDProb, Prob};
use bv::BitVec;
use datavzrd::render_report;
use serde_yaml::Value;

use derefable::Derefable;

use derive_deref::DerefMut;

use ordered_float::NotNan;

use rust_htslib::bcf::{
    self,
    record::GenotypeAllele::{Phased, Unphased},
    Read,
};

use good_lp::IntoAffineExpression;
use good_lp::*;
use good_lp::{variable, Expression};

use petgraph::dot::{Config, Dot};
use petgraph::graph::{Graph, NodeIndex};

use serde::Serialize;
use serde_json::json;
use std::collections::{BTreeMap, HashMap};

use derive_new::new;
use std::collections::BTreeSet;
use std::fs;
use std::io::Write;
use std::str::FromStr;
use std::{path::PathBuf, str};

#[derive(Derefable, Debug, Copy, Clone, PartialEq, Eq, Hash, Ord, PartialOrd, Serialize)]
pub struct VariantID(#[deref] pub i32);

#[derive(Derefable, Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize)]
pub struct Haplotype(#[deref] pub String);

#[derive(Debug, Clone, Derefable)]
pub struct AlleleFreqDist(#[deref] BTreeMap<AlleleFreq, LogProb>);

impl AlleleFreqDist {
    pub fn vaf_query(&self, vaf: &AlleleFreq) -> Option<LogProb> {
        if self.contains_key(vaf) {
            return Some(*self.get(vaf).unwrap());
        } else {
            let (x_0, y_0) = self.range::<NotNan<f64>, _>(..vaf).next_back().unwrap();
            let (x_1, y_1) = self.range::<NotNan<f64>, _>(vaf..).next().unwrap();
            // If y_0 == y_1, what to do?
            // if y_0 == y_1 {
            //     return Some(*y_0);
            // }
            //calculate the proportion (vaf - x_0) / (x_1 - x_0) in linear space
            let fraction = (*vaf - *x_0) / (*x_1 - *x_0);

            //ensure we interpolate in the right direction
            let (smaller, larger) = if y_0 < y_1 { (y_0, y_1) } else { (y_1, y_0) };

            //take log difference (always positive) and subtract from each other
            let log_diff = larger.ln_sub_exp(*smaller).exp();

            //interpolated probability in linear space again
            let interpolated_prob = smaller.exp() + fraction * log_diff;

            //convert back to LogProb
            Some(LogProb::from(Prob(interpolated_prob)))
        }
    }
}

#[derive(Derefable, Debug, Clone)]
pub struct CandidateMatrix(#[deref] BTreeMap<VariantID, (BitVec, BitVec)>);

impl CandidateMatrix {
    pub fn new(
        haplotype_variants: &HaplotypeVariants,
        // haplotypes: &Vec<Haplotype>,
    ) -> Result<Self> {
        let mut candidate_matrix = BTreeMap::new();
        haplotype_variants.iter().for_each(|(variant_id, bmap)| {
            let mut haplotype_variants_gt = BitVec::new();
            let mut haplotype_variants_c = BitVec::new();
            bmap.iter().for_each(|(_haplotype, (gt, c))| {
                haplotype_variants_c.push(*c);
                haplotype_variants_gt.push(*gt);
            });
            candidate_matrix.insert(*variant_id, (haplotype_variants_gt, haplotype_variants_c));
        });
        Ok(CandidateMatrix(candidate_matrix))
    }
    pub fn find_identical_haplotypes(
        &self,
        haplotypes: Vec<Haplotype>,
    ) -> (
        BTreeMap<Haplotype, Vec<Haplotype>>,
        BTreeMap<Haplotype, Vec<Haplotype>>,
    ) {
        let mut signature_map: BTreeMap<Vec<bool>, Vec<Haplotype>> = BTreeMap::new();

        // Iterate over each haplotype to generate its signature
        for (index, haplotype) in haplotypes.iter().enumerate() {
            let mut signature: Vec<bool> = Vec::new();

            // Create the signature by extracting bits from each variant's BitVec
            for (_variant_id, (bitvec, _)) in &self.0 {
                let presence = bitvec.get(index as u64);
                signature.push(presence);
            }

            // Group haplotypes with the same signature
            signature_map
                .entry(signature)
                .or_insert_with(Vec::new)
                .push(haplotype.clone());
        }

        // Create a BTreeMap where each haplotype from groups has its own entry
        let mut haplotype_map: BTreeMap<Haplotype, Vec<Haplotype>> = BTreeMap::new();

        //and create a BTreeMap where every group has only one representative key
        let mut haplotype_map_representative: BTreeMap<Haplotype, Vec<Haplotype>> = BTreeMap::new();

        for (_signature, group) in signature_map {
            let representative = group[0].clone(); // First haplotype as representative
            haplotype_map_representative.insert(representative, group.clone());

            for haplotype in &group {
                haplotype_map.insert(haplotype.clone(), group.clone());
            }
        }

        (haplotype_map_representative, haplotype_map)
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
pub struct VariantCalls(#[deref] BTreeMap<VariantID, (String, f32, AlleleFreqDist, i32)>); //The place of f32 is maximum a posteriori estimate of AF.

impl VariantCalls {
    pub fn new(variant_calls: &mut bcf::Reader) -> Result<Self> {
        let mut calls = BTreeMap::new();
        let header = variant_calls.header().clone();
        for record_result in variant_calls.records() {
            let mut record = record_result?;
            record.unpack();
            let prob_absent = record.info(b"PROB_ABSENT").float().unwrap().unwrap()[0];
            let _prob_absent_prob = Prob::from(PHREDProb(prob_absent.into()));
            let afd_utf = record.format(b"AFD").string()?;
            let afd = std::str::from_utf8(afd_utf[0]).unwrap();
            let read_depth = record.format(b"DP").integer().unwrap();
            let read_depth_int = read_depth[0].get(0).unwrap();
            //include all variants without making a difference for their depth of coverage.
            // if read_depth[0] != &[0]
            // // && (&prob_absent_prob <= &Prob(0.05) || &prob_absent_prob >= &Prob(0.95))
            // {
            //because some afd strings are just "." and that throws an error while splitting below.
            let variant_id: i32 = String::from_utf8(record.id())?.parse().unwrap();
            let af = (&*record.format(b"AF").float().unwrap()[0]).to_vec()[0];
            let mut vaf_density = BTreeMap::new();
            for pair in afd.split(',') {
                if let Some((vaf, density)) = pair.split_once("=") {
                    let (vaf, density): (AlleleFreq, f64) =
                        (vaf.parse().unwrap(), density.parse().unwrap());
                    vaf_density.insert(vaf, LogProb::from(PHREDProb(density)));
                }
            }
            //exxtract reference (REF), alternative (ALT), and position (POS)
            let chr_name =
                std::str::from_utf8(header.rid2name(record.rid().unwrap()).unwrap()).unwrap();
            let pos = record.pos() + 1; // 1-based indexing
            let ref_base = std::str::from_utf8(record.alleles()[0]).unwrap();
            let alt_base = std::str::from_utf8(record.alleles()[1]).unwrap(); // Only first ALT allele

            //construct the variant change string (e.g., "C1973T")
            let variant_change = format!("{}:{}{}{}", chr_name, ref_base, pos, alt_base);

            calls.insert(
                VariantID(variant_id),
                (
                    variant_change,
                    af,
                    AlleleFreqDist(vaf_density),
                    *read_depth_int,
                ),
            );
            // }
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
    pub fn without_zero_dp(&self) -> VariantCalls {
        // Filter out entries with 0 in the last element of the tuple
        let filtered_map: BTreeMap<_, _> = self
            .0
            .iter()
            .filter(|(_, &(_, _, _, dp))| dp != 0) // Keep only entries where dp != 0
            .map(|(k, v)| (k.clone(), v.clone())) // Clone key-value pairs to create a new map
            .collect();

        VariantCalls(filtered_map) // Return a new VariantCalls struct
    }
    pub fn check_variant_threshold(
        &self,
        haplotype_variants: &HaplotypeVariants,
        threshold_considered_variants: f64,
    ) -> Result<bool> {
        let variants_haplotype_variants: Vec<_> = haplotype_variants.keys().cloned().collect();
        let variants_haplotype_calls: Vec<_> = self.keys().cloned().collect();

        let rateof_evaluated_variants: f64 =
            variants_haplotype_calls.len() as f64 / variants_haplotype_variants.len() as f64;
        Ok(rateof_evaluated_variants > threshold_considered_variants)
    }
}

#[derive(Derefable, Debug, Clone, PartialEq, Eq, PartialOrd, DerefMut)]
pub struct HaplotypeVariants(#[deref] pub BTreeMap<VariantID, BTreeMap<Haplotype, (bool, bool)>>);

#[derive(Derefable, Debug, Clone, new)]
pub struct HaplotypeGraph {
    #[deref]
    graph: Graph<(Haplotype, Haplotype), i32, petgraph::Undirected>,
    node_indices: HashMap<(Haplotype, Haplotype), NodeIndex>,
}

impl HaplotypeGraph {
    pub fn get_node_index(&self, query: &(Haplotype, Haplotype)) -> Option<NodeIndex> {
        self.node_indices.get(query).cloned()
    }
}

#[derive(Derefable, Debug, Clone)]
pub struct DistanceMatrix(#[deref] pub BTreeMap<(Haplotype, Haplotype), usize>);

impl HaplotypeVariants {
    pub fn new(haplotype_variants: &mut bcf::Reader) -> Result<Self> {
        let mut variant_records = BTreeMap::new();
        for record_result in haplotype_variants.records() {
            let record = record_result?;
            let variant_id: VariantID = VariantID(String::from_utf8(record.id())?.parse().unwrap());
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
                        matrices.insert(haplotype.clone(), (true, loci[index] == &[1]));
                    } else {
                        matrices.insert(haplotype.clone(), (false, loci[index] == &[1]));
                    }
                }
            }
            variant_records.insert(variant_id, matrices);
        }
        Ok(HaplotypeVariants(variant_records))
    }

    pub fn filter_for_variants(&self, variant_ids: &Vec<VariantID>) -> Result<HaplotypeVariants> {
        let mut filtered_haplotype_variants: BTreeMap<
            VariantID,
            BTreeMap<Haplotype, (bool, bool)>,
        > = BTreeMap::new();
        for (variant, haplotype_map) in self.iter() {
            if variant_ids.contains(&variant) {
                filtered_haplotype_variants.insert(variant.clone(), haplotype_map.clone());
            }
        }
        Ok(HaplotypeVariants(filtered_haplotype_variants))
    }

    pub fn filter_for_haplotypes(&self, haplotypes: &Vec<Haplotype>) -> Result<Self> {
        let mut new_haplotype_variants: BTreeMap<VariantID, BTreeMap<Haplotype, (bool, bool)>> =
            BTreeMap::new();
        for (variant, matrix_map) in self.iter() {
            let mut new_matrix_map = BTreeMap::new();
            for (haplotype_m, (variant_status, coverage_status)) in matrix_map {
                if haplotypes.contains(&haplotype_m) {
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
        let candidate_matrix_values: Vec<(BitVec, BitVec)> = CandidateMatrix::new(self)
            .unwrap()
            .values()
            .cloned()
            .collect();
        let mut common_variants = Vec::new();
        for ((_genotype_matrix, coverage_matrix), (variant, (_, _, _, _))) in
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

    pub fn find_equivalence_classes_with_graph(
        &self,
        application: &str,
        threshold: usize, //an edge in the graph representation for the equivalence classes is drawn if and only if the distance in terms of variants is smaller than a given threshold and the two nodes belong to the same group
        output_graph: &PathBuf,
    ) -> Result<HaplotypeGraph> {
        //create file path  for the graph
        let mut parent = output_graph.clone();
        parent.pop();
        let mut file = fs::File::create(parent.join("graph.dot")).unwrap();

        let mut equivalence_classes: BTreeMap<Haplotype, Vec<VariantID>> = BTreeMap::new();

        //initialize equivalence classes with haplotypes (important for haplotypes that have no variant e.g. A*03:01)
        self.iter().for_each(|(_, haplotype_map)| {
            haplotype_map.iter().for_each(|(haplotype, _)| {
                equivalence_classes.insert(haplotype.clone(), vec![]);
            })
        });

        for (variant, haplotype_map) in self.iter() {
            for (haplotype, (variant_in_gt, _variant_in_c)) in haplotype_map.iter() {
                match variant_in_gt {
                    true => {
                        let mut variants_in_haplotype = equivalence_classes[&haplotype].clone();
                        variants_in_haplotype.push(variant.clone());
                        equivalence_classes.insert(haplotype.clone(), variants_in_haplotype);
                    }
                    _ => (),
                }
            }
        }

        let mut deps = Graph::new_undirected();

        for (haplotype, variants) in equivalence_classes.iter() {
            //initialize variables
            let mut splitted_1 = vec![];
            let mut haplotype_group = Haplotype(String::from(""));
            if &application == &"hla" {
                splitted_1 = haplotype.split(':').collect::<Vec<&str>>();
                haplotype_group = Haplotype(splitted_1[0].to_owned() + &":" + splitted_1[1]);
            } else if &application == &"virus" {
                haplotype_group = haplotype.clone();
            }

            //add new node
            deps.add_node((haplotype.clone(), haplotype_group.clone()));

            //find index of the main node. this is necessary for the check and drawing an edge later
            let index = deps
                .node_indices()
                .find(|i| deps[*i] == (haplotype.clone(), haplotype_group.clone()))
                .unwrap();

            //loop over the graph to draw edges if conditions are met
            for idx in 0..deps.node_count() {
                // find the node and its variants at index
                let node_at_index = &deps[NodeIndex::new(idx)];
                let variants_of_node_at_index: &Vec<VariantID> =
                    &equivalence_classes[&node_at_index.0];

                let mut splitted_2 = vec![];
                let mut haplotype_group_at_index = Haplotype(String::from(""));
                if &application == &"hla" {
                    splitted_2 = node_at_index.0.split(':').collect::<Vec<&str>>();
                    haplotype_group_at_index =
                        Haplotype(splitted_2[0].to_owned() + &":" + splitted_2[1]);
                } else if &application == &"virus" {
                    haplotype_group_at_index = node_at_index.0.clone();
                }

                //find the differences between the main node and the neighbor node
                let mut difference: Vec<&VariantID> = vec![];
                for variant in variants_of_node_at_index.iter() {
                    if !variants.contains(&variant) {
                        difference.push(variant);
                    }
                }
                //for hla: draw an edge only if the difference is less than the threshold, haplotype group is the same
                // and it's not the same node
                //for virus: draw an edge only if the difference is less than the threshold and it's not the same node

                if &application == &"hla"
                    && (difference.len() < threshold)
                    && (haplotype_group == haplotype_group_at_index)
                    && (index != NodeIndex::new(idx))
                {
                    deps.add_edge(index, NodeIndex::new(idx), 1);
                } else if &application == &"virus"
                    && (difference.len() < threshold)
                    && (index != NodeIndex::new(idx))
                {
                    deps.add_edge(index, NodeIndex::new(idx), 1);
                }
            }
        }

        //find node indices
        let mut node_indices = HashMap::new();

        //populate node_indices by iterating through all nodes
        for node_index in deps.node_indices() {
            let data = deps[node_index].clone(); // Get the data of the node
            node_indices.insert(data, node_index);
        }

        //write graph to dot file path
        let output = format!("{:?}", Dot::with_config(&deps, &[Config::EdgeNoLabel]));
        file.write_all(&output.as_bytes())
            .expect("could not write file");
        Ok(HaplotypeGraph {
            graph: deps,
            node_indices: node_indices,
        })
    }

    pub fn find_equivalence_classes_hamming_distance(
        &self,
        _application: &str,
    ) -> Result<DistanceMatrix> {
        //initialize distances as a HashMap
        let mut distances: BTreeMap<(Haplotype, Haplotype), usize> = BTreeMap::new();

        //compute hamming distances between each pair of haplotypes
        let haplotype_keys: Vec<Haplotype> = self
            .values()
            .cloned()
            .collect::<Vec<BTreeMap<Haplotype, (bool, bool)>>>()[0]
            .keys()
            .cloned()
            .collect();
        for i in 0..haplotype_keys.len() {
            for j in i + 1..haplotype_keys.len() {
                let mut distance = 0;
                let hap1 = &haplotype_keys[i];
                let hap2 = &haplotype_keys[j];
                for (_variant, haplotype_map) in self.iter() {
                    let gt1 = &haplotype_map[&hap1].0;
                    let gt2 = &haplotype_map[&hap2].0;
                    if *gt1 != *gt2 {
                        distance += 1;
                    }
                }
                //insert the distance into the distances HashMap
                distances.insert((hap1.clone(), hap2.clone()), distance);
            }
        }
        //print the distances
        for ((hap1, hap2), distance) in &distances {
            println!("Distance between {:?} and {:?}: {}", hap1, hap2, distance);
        }
        Ok(DistanceMatrix(distances))
    }
}

#[derive(Serialize, Debug)]
pub(crate) struct DatasetVariants {
    variant_change: String,
    vaf: f32,
}
#[derive(Serialize, Debug)]
pub(crate) struct DatasetHaplotypeVariants {
    variant_change: String,
    haplotype: String,
}
#[derive(Serialize, Debug)]
pub(crate) struct DatasetHaplotypeFractions {
    haplotype: String,
    fraction: AlleleFreq,
}

#[derive(Serialize, Debug)]
pub(crate) struct DatasetAfd {
    variant_change: String,
    allele_freq: AlleleFreq,
    probability: f64,
}

pub fn plot_prediction(
    outdir: &PathBuf,
    solution: &str,
    candidate_matrix_values: &Vec<(BitVec, BitVec)>,
    haplotypes: &Vec<Haplotype>,
    variant_calls: &VariantCalls,
    best_variables: &Vec<f64>,
) -> Result<()> {
    let mut parent = outdir.clone();
    parent.pop();
    fs::create_dir_all(&parent)?;

    let mut file_name = "".to_string();
    let mut json = include_str!("../../../templates/final_prediction.json");
    if &solution == &"lp" {
        json = include_str!("../../../templates/lp_prediction.json");
    }
    let mut blueprint: serde_json::Value = serde_json::from_str(json).unwrap();
    let mut plot_data_variants = Vec::new();
    let mut plot_data_haplotype_variants = Vec::new();
    let mut plot_data_haplotype_fractions = Vec::new();
    let mut plot_data_covered_variants = Vec::new();
    let mut plot_data_dataset_afd = Vec::new();

    if &solution == &"lp" {
        //write tsv table for datavzrd view
        let mut variant_records = Vec::new();
        let lp_solution_path = parent.join("lp_solution.tsv");
        let mut wtr_lp = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .quote_style(csv::QuoteStyle::Never)
            .from_path(&lp_solution_path)?;

        //filter haplotypes based on any variant they contain has nonzero AF
        let mut contributing_haplotypes = std::collections::HashSet::new();
        for ((genotype_matrix, coverage_matrix), (_, (var_change, af, _, _dp))) in
            candidate_matrix_values.iter().zip(variant_calls.iter())
        {
            if *af > 0.0 {
                for (i, haplotype) in haplotypes.iter().enumerate() {
                    let has_variant = genotype_matrix[i as u64];
                    if has_variant {

                        contributing_haplotypes.insert(haplotype.to_string()); // NEW
                    }
                }
            }
        }
        dbg!(&contributing_haplotypes.len());

        // print haplotypes that will not be displayed
        // for haplotype in haplotypes.iter() {
        //     if !contributing_haplotypes.contains(&haplotype.to_string()) {
        //         dbg!("haplotype not in contributing_haplotypes: {}", &haplotype);
        //     }
        // }

        //filter haplotypes and variables based on contributing haplotypes
        let filtered: Vec<_> = haplotypes.iter()
            .zip(best_variables.iter())
            .filter(|(h, _)| contributing_haplotypes.contains(&h.to_string())) 
            .collect();

        let (filtered_haplotypes, filtered_best_variables): (Vec<_>, Vec<_>) = filtered.into_iter().unzip(); 

        //sort haplotypes descending by lp values
        let mut haplotype_with_weights: Vec<(&Haplotype, &f64)> =
        filtered_haplotypes
            .iter()
            .zip(filtered_best_variables.iter())
            .map(|(&h, &v)| (h, v))
            .collect();
        haplotype_with_weights.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap_or(std::cmp::Ordering::Equal));

        let sorted_haplotypes: Vec<String> = haplotype_with_weights.iter().map(|(h, _)| (*h).to_string()).collect();
        let sorted_indices: Vec<usize> = haplotype_with_weights
            .iter()
            .map(|(h, _)| filtered_haplotypes.iter().position(|x| *x == *h).unwrap())
            .collect();

        let filtered_candidate_matrix_values: Vec<(bv::BitVec, bv::BitVec)> = candidate_matrix_values
            .iter()
            .map(|(genotype_matrix, coverage_matrix)| {
                let mut filtered_genotype = bv::BitVec::new();
                let mut filtered_coverage = bv::BitVec::new();

                for (i, haplotype) in haplotypes.iter().enumerate() {
                    if contributing_haplotypes.contains(&haplotype.to_string()) {
                        filtered_genotype.push(genotype_matrix[i as u64]);
                        filtered_coverage.push(coverage_matrix[i as u64]);
                    }
                }

                (filtered_genotype, filtered_coverage)
            })
            .collect();

        let mut headers: Vec<_> = vec![
            "sum_of_fractions".to_string(),
            "variant".to_string(),
            "vaf".to_string(),
            "prediction_error".to_string(),
        ];
        //insert sorted haplotype names to header
        headers.extend(sorted_haplotypes.clone());
        wtr_lp.write_record(&headers)?;

        for ((genotype_matrix, coverage_matrix), (variant_id, (var_change, af, _, _dp))) in
        filtered_candidate_matrix_values.iter().zip(variant_calls.iter())
        {
            let mut counter = 0;
            for (i, _variable) in filtered_best_variables.iter().enumerate() {
                if coverage_matrix[i as u64] {
                    counter += 1;
                }
            }
            if counter == filtered_best_variables.len() {
                //initialize the dict for storing sum_of_fractions (required for the table)
                let mut haplotype_has_variant = Vec::new();
                let mut sum_of_fractions_vec = Vec::new();

                for (i, (variable, haplotype)) in
                filtered_best_variables.iter().zip(filtered_haplotypes.iter()).enumerate()
                {
                    //record presence/absence for the variant
                    let has_variant = genotype_matrix[i as u64];
                    // haplotype_has_variant.push(has_variant.to_string());
                    haplotype_has_variant.push(if has_variant { "1" } else { "0" }.to_string());
                    
                    if has_variant {
                        plot_data_haplotype_fractions.push(DatasetHaplotypeFractions {
                            haplotype: haplotype.to_string(),
                            fraction: NotNan::new(**variable).unwrap(),
                        });
                        plot_data_haplotype_variants.push(DatasetHaplotypeVariants {
                            variant_change: var_change.to_string(),
                            haplotype: haplotype.to_string(),
                        });
                        plot_data_variants.push(DatasetVariants {
                            variant_change: var_change.to_string(),
                            vaf: af.clone(),
                        });

                        //fill in the dict required for the table
                        let sum_of_fractions_map = json!({
                            "haplotype": haplotype.to_string(),
                            "vaf": af,
                            "fraction": variable
                        });
                        sum_of_fractions_vec.push(sum_of_fractions_map);
                    }
                }
                //push other members and write to table if at least one haplotype contains the variant.
                if !sum_of_fractions_vec.is_empty() {
                    let json_array_string = serde_json::to_string(&sum_of_fractions_vec)?;
                    let mut final_record = Vec::new();
                    final_record.push(json_array_string);
                    final_record.push(var_change.to_string());

                    //push vaf
                    final_record.push(af.to_string());

                    //sum all fractions
                    let total_fraction: f64 = sum_of_fractions_vec
                        .iter()
                        .filter_map(|entry| entry.get("fraction"))
                        .filter_map(|val| val.as_f64())
                        .sum();

                    //find prediction error and push
                    let difference = (af - total_fraction as f32).abs();
                    final_record.push(difference.to_string());

                    //push haplotype presence information
                    let sorted_variant_flags: Vec<String> = sorted_indices
                        .iter()
                        .map(|&i| haplotype_has_variant[i].clone())
                        .collect();
                    final_record.extend(sorted_variant_flags);
                    variant_records.push(final_record);
                }
            }
        }

        //then write each row to tsv
        for record in variant_records {
            wtr_lp.write_record(&record)?;
        }
        // render datavzrd report
        let contents = include_str!("../../../templates/datavzrd_config.yaml");

        //parse YAML into a serde_yaml::Value
        let mut config_yaml: Value = serde_yaml::from_str(&contents)?;

        //navigate through the yaml and change path
        let lp_solution_path_str = lp_solution_path.into_os_string().into_string().unwrap();
        dbg!(&lp_solution_path_str);
        if let Some(datasets) = config_yaml.get_mut("datasets") {
            if let Some(solutions) = datasets.get_mut("solutions") {
                if let Some(path_field) = solutions.get_mut("path") {
                    *path_field = Value::String(lp_solution_path_str);
                } else {
                    eprintln!("Warning: 'path' field not found under 'datasets.solutions'");
                }
            } else {
                eprintln!("Warning: 'solutions' not found under 'datasets'");
            }
        } else {
            eprintln!("Warning: 'datasets' not found in config");
        }
        //write updated config back to a new file (or overwrite the original)
        let filled_config_yaml = parent.join("datavzrd_config_filled.yaml");
        let yaml_output = PathBuf::from(&filled_config_yaml);
        fs::write(&yaml_output, serde_yaml::to_string(&config_yaml)?)?;

        //then create datavzrd report
        let datavzrd_output = parent.join("datavzrd_report");
        render_report(&filled_config_yaml, &datavzrd_output, "", false, true)?;

        file_name.push_str("lp_solution.json");
        wtr_lp.flush()?;
    } else if &solution == &"final" {
        candidate_matrix_values
            .iter()
            .zip(variant_calls.iter())
            .for_each(
                |((genotypes, covered), (variant_id, (var_change, af, afd, dp)))| {
                    let mut b_check = false;
                    let mut b_fraction = 0.0;
                    best_variables
                        .iter()
                        .zip(haplotypes.iter())
                        .enumerate()
                        .for_each(|(i, (fraction, haplotype))| {
                            //create an exception for B to add it both to haplotype variants and fractions
                            if haplotype == &Haplotype("B".to_string()) {
                                b_fraction = *fraction;
                            }
                            if genotypes[i as u64] && covered[i as u64] {
                                if *dp != 0 {
                                    plot_data_variants.push(DatasetVariants {
                                        variant_change: var_change.to_string(),
                                        vaf: *af,
                                    });
                                }
                                plot_data_haplotype_fractions.push(DatasetHaplotypeFractions {
                                    haplotype: haplotype.to_string(),
                                    fraction: NotNan::new(*fraction).unwrap(),
                                });
                                plot_data_haplotype_variants.push(DatasetHaplotypeVariants {
                                    variant_change: var_change.to_string(),
                                    haplotype: haplotype.to_string(),
                                });
                                for (j, haplotype) in haplotypes.iter().enumerate() {
                                    if covered[j as u64] {
                                        plot_data_covered_variants.push(DatasetHaplotypeVariants {
                                            variant_change: var_change.to_string(),
                                            haplotype: haplotype.to_string(),
                                        });
                                    }
                                }
                                //also add the heatmap for afd below the covered panels
                                for (allele_freq, prob) in afd.iter() {
                                    plot_data_dataset_afd.push(DatasetAfd {
                                        variant_change: var_change.to_string(),
                                        allele_freq: *allele_freq,
                                        probability: f64::from(*prob),
                                    })
                                }
                                b_check = true;
                            } else {
                                if *af > 0.0 {
                                    if *dp != 0 {
                                        plot_data_variants.push(DatasetVariants {
                                            variant_change: var_change.to_string(),
                                            vaf: *af,
                                        });
                                    }
                                    for (j, haplotype) in haplotypes.iter().enumerate() {
                                        if covered[j as u64] {
                                            plot_data_covered_variants.push(
                                                DatasetHaplotypeVariants {
                                                    variant_change: var_change.to_string(),
                                                    haplotype: haplotype.to_string(),
                                                },
                                            );
                                        }
                                    }
                                    //also add the heatmap for afd below the covered panels
                                    for (allele_freq, prob) in afd.iter() {
                                        plot_data_dataset_afd.push(DatasetAfd {
                                            variant_change: var_change.to_string(),
                                            allele_freq: *allele_freq,
                                            probability: f64::from(*prob),
                                        })
                                    }
                                }
                            }
                        });
                    if b_check && haplotypes.contains(&Haplotype("B".to_string())) {
                        //todo: refactor this as this is not ideal to check for every sample and every variant.
                        plot_data_haplotype_variants.push(DatasetHaplotypeVariants {
                            variant_change: var_change.to_string(),
                            haplotype: "B".to_string(),
                        });
                        plot_data_haplotype_fractions.push(DatasetHaplotypeFractions {
                            haplotype: "B".to_string(),
                            fraction: NotNan::new(b_fraction).unwrap(),
                        });
                    }
                },
            );
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
    extend_haplotypes: bool,
    num_variant_distance: i64,
    num_constraint_haplotypes: i32,
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

    // introduce auxiliary binary variables for sparsity (these will act as switches for real variables)
    let binary_vars: Vec<Variable> =
        problem.add_vector(variable().integer().min(0.0).max(1.0), variables.len());

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

    // add constraints linking binary variables and the original variables
    let m_constant = 100;
    for (var, bin_var) in variables.iter().zip(binary_vars.iter()) {
        model = model.with(constraint!(*bin_var >= *var));
        model = model.with(constraint!(*bin_var <= m_constant * *var));
    }

    // add the sparsity constraint: sum of binary variables <= 6
    let mut binary_sum = Expression::from_other_affine(0);
    for bin_var in binary_vars.iter() {
        binary_sum += Expression::from_other_affine(bin_var);
    }
    model = model.with(constraint!(binary_sum == num_constraint_haplotypes));

    //add the constraints to the model
    for (c, t_var) in constraints.iter().zip(t_vars.iter()) {
        model = model.with(constraint!(t_var >= c.clone()));
        model = model.with(constraint!(t_var >= -c.clone()));
    }

    //solve the problem with the default solver, i.e. coin_cbc
    let solution = model.solve().unwrap();

    let mut best_variables = Vec::new();
    //finally, print the variables and the sum
    let mut lp_haplotypes = BTreeMap::new();
    for (i, (var, haplotype)) in variables.iter().zip(haplotypes.iter()).enumerate() {
        println!("v{}, {}={}", i, haplotype.to_string(), solution.value(*var));
        best_variables.push(solution.value(var.clone()).clone());
        if solution.value(*var) > lp_cutoff {
            //the speed of fraction exploration is managable in case of diploid priors
            lp_haplotypes.insert(haplotype.clone(), solution.value(*var).clone());
        }
    }
    // println!("sum = {}", solution.eval(sum_tvars));
    //plot the best result
    let candidate_matrix_values: Vec<(BitVec, BitVec)> =
        candidate_matrix.values().cloned().collect();
    plot_prediction(
        outdir,
        &"lp",
        &candidate_matrix_values,
        &haplotypes,
        &variant_calls,
        &best_variables,
    )?;
    let lp_haplotypes_keys: Vec<_> = lp_haplotypes.keys().cloned().collect();
    //extend haplotypes found by linear program, add haplotypes that have the same variants to the final list.
    //then sort by hamming distance, take the closest x additional alleles according to 'num_variant_distance'.
    //this is done using haplotype_dict, hence ONLY the nonzero DP variants that are covered by all haplotypes (C:1) are included in this extension.
    let mut extended_haplotypes_bset: BTreeSet<_> = lp_haplotypes_keys.iter().cloned().collect();
    if extend_haplotypes {
        lp_haplotypes.iter().for_each(|(f_haplotype, _)| {
            let variants = haplotype_dict.get(&f_haplotype).unwrap().clone();
            let variant_set: BTreeSet<_> = variants.iter().collect();
            haplotype_dict
                .iter()
                .for_each(|(haplotype, haplotype_variants)| {
                    let haplotype_set: BTreeSet<_> = haplotype_variants.iter().collect();
                    let difference: Vec<_> =
                        variant_set.symmetric_difference(&haplotype_set).collect();
                    //0 distance is excluded because lp gets only the representative, nonidentical haplotypes as input.
                    if (difference.len() as i64 <= num_variant_distance
                        && difference.len() as i64 != 0)
                        && ((variants.len() as i64 - haplotype_variants.len() as i64).abs()
                            <= num_variant_distance)
                    //the last check is necessary to make sure the difference on haplotypes is also correlated with the lengths of vectors of variants.
                    {
                        extended_haplotypes_bset.insert(haplotype.clone());
                    }
                });
        });
        let extended_haplotypes = extended_haplotypes_bset.into_iter().collect();
        dbg!(&extended_haplotypes);
        Ok(extended_haplotypes)
    } else {
        Ok(lp_haplotypes_keys)
    }
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
        .map(|(_, (_, _, afd, _))| afd.clone())
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
                        if genotypes[i as u64] && covered[i as u64] {
                            vaf_sum += *fraction;
                            counter += 1;
                        } else if genotypes[i as u64] == false && covered[i as u64] == false {
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
    let best_odds: f64 = 1.00;
    let format_f64 = |number: f64, records: &mut Vec<String>| {
        if number <= 0.01 {
            records.push(format!("{:+.2e}", number))
        } else {
            records.push(format!("{:.2}", number))
        }
    };
    format_f64(best_density.exp(), &mut records);
    records.push(format!("{:.2}", best_odds));
    let format_freqs = |frequency: NotNan<f64>, records: &mut Vec<String>| {
        if frequency <= NotNan::new(0.01).unwrap() {
            records.push(format!("{:+.4e}", NotNan::into_inner(frequency)))
        } else {
            records.push(format!("{:.4}", frequency))
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
    let candidate_matrix_values: Vec<(BitVec, BitVec)> =
        candidate_matrix.values().cloned().collect();
    //collect haplotype-to-variants information
    let mut haplotype_dict: HashMap<Haplotype, Vec<VariantID>> =
        haplotypes.iter().map(|h| (h.clone(), vec![])).collect();
    //variant-wise iteration
    let mut expr = Expression::from_other_affine(0.); // A constant expression
    for ((genotype_matrix, coverage_matrix), (variant, (_, af, _, _))) in
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
                if genotype_matrix[i as u64] {
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

pub fn get_event_posteriors(
    haplotype_variants: &HaplotypeVariants,
    variant_calls: VariantCalls,
    application: &str,
    prior: &String,
    outfile: &PathBuf,
    extend_haplotypes: bool,
    num_extend_haplotypes: i64,
    num_constraint_haplotypes: i32,
    lp_cutoff: f64,
    enable_equivalence_class_constraint: bool,
    threshold_equivalence_class: Option<usize>,
) -> Result<(Vec<(HaplotypeFractions, LogProb)>, Vec<Haplotype>, Data)> {
    //FIRST, perform linear program using only nonzero DP variants
    //filter variant calls and haplotype variants
    let filtered_calls = variant_calls.without_zero_dp();
    let nonzero_dp_variants: Vec<VariantID> = filtered_calls.keys().cloned().collect();
    let var_filt_haplotype_variants: HaplotypeVariants =
        haplotype_variants.filter_for_variants(&nonzero_dp_variants)?;

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
    // Print the result
    // for (representative, group) in &identical_haplotypes_map {
    //     println!("Representative: {:?}, Group: {:?}", representative, group);
    // }
    let representatives = identical_haplotypes_map_rep.keys().cloned().collect();
    let repr_haplotype_variants =
        var_filt_haplotype_variants.filter_for_haplotypes(&representatives)?;
    let repr_candidate_matrix = CandidateMatrix::new(&repr_haplotype_variants).unwrap();

    //employ the linear program
    let lp_haplotypes = linear_program(
        &outfile,
        &repr_candidate_matrix,
        &representatives,
        &filtered_calls,
        lp_cutoff,
        extend_haplotypes,
        num_extend_haplotypes, //extension functionality for virus case is not recommended for now as it will lead to performance problems.
        num_constraint_haplotypes,
    )?;

    //SECOND, model evaluation using ALL variants but only the LP- selected haplotypes
    //prepare inputs of model evaluation
    let hap_filt_haplotype_variants = haplotype_variants.filter_for_haplotypes(&lp_haplotypes)?;
    let model_candidate_matrix = CandidateMatrix::new(&hap_filt_haplotype_variants)?;

    //compute model
    let prior = PriorTypes::from_str(&prior).unwrap();
    let model = Model::new(
        Likelihood::new(),
        Prior::new(prior.clone()),
        Posterior::new(),
    );
    let data = Data::new(model_candidate_matrix.clone(), variant_calls.clone());

    let mut eq_graph: Option<_> = Some(HaplotypeGraph::new(Graph::default(), HashMap::new()));
    let mut application_name = "none";

    if &application == &"hla" {
        //equivalence graph based optimization is at the developmental phase.
        eq_graph = Some(
            hap_filt_haplotype_variants
                .find_equivalence_classes_with_graph(
                    "hla",
                    threshold_equivalence_class.unwrap(),
                    &outfile,
                )
                .unwrap(),
        );
        application_name = &"hla";
    } else if &application == &"virus" {
        eq_graph = None;
        application_name = &"virus";
    }
    let computed_model = model.compute_from_marginal(
        &Marginal::new(
            lp_haplotypes.len(),
            lp_haplotypes.clone(),
            prior.clone(),
            eq_graph,
            enable_equivalence_class_constraint,
            application_name.to_string(),
        ),
        &data,
    );

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
        &prior,
        &identical_haplotypes_map,
    )?;
    Ok((new_event_posteriors, all_haplotypes, data))
}

fn generate_combinations(
    expanded_fractions: &Vec<AlleleFreq>,
    haplotype_indices: &BTreeMap<Haplotype, usize>,
    prior: &PriorTypes,
    identical_haplotypes_map: &BTreeMap<Haplotype, Vec<Haplotype>>,
) -> Vec<Vec<AlleleFreq>> {
    let mut result = vec![expanded_fractions.clone()];

    for (cur_haplotype, &idx_current) in haplotype_indices.iter() {
        let fraction = expanded_fractions[idx_current];

        if fraction > NotNan::new(0.0).unwrap() {
            if let Some(identical_haplotypes) = identical_haplotypes_map.get(cur_haplotype) {
                let mut new_combinations = Vec::new();
                let mut previous_haplotypes = Vec::new(); // Store previous haplotypes in the group

                for existing_row in &result {
                    for ident_h in identical_haplotypes {
                        if ident_h != cur_haplotype {
                            let idx_ident_h = haplotype_indices.get(ident_h).unwrap();

                            //standard case: transfer full fraction
                            let mut alt_row = existing_row.clone();
                            alt_row[*idx_ident_h] = fraction;
                            alt_row[idx_current] = NotNan::new(0.0).unwrap();
                            new_combinations.push(alt_row);

                            //diploid-specific case: If fraction is 1.0, create a row with 0.5 fractions
                            if let PriorTypes::Diploid = prior {
                                if fraction == NotNan::new(1.0).unwrap() {
                                    let mut diploid_row = existing_row.clone();
                                    diploid_row[*idx_ident_h] = NotNan::new(0.5).unwrap();
                                    diploid_row[idx_current] = NotNan::new(0.5).unwrap();
                                    new_combinations.push(diploid_row);

                                    //handle within-group pairs in diploid case
                                    for prev_hap in &previous_haplotypes {
                                        let idx_prev = haplotype_indices.get(prev_hap).unwrap();

                                        if idx_ident_h != idx_prev {
                                            let mut pair_row = existing_row.clone();
                                            pair_row[*idx_ident_h] = NotNan::new(0.5).unwrap();
                                            pair_row[*idx_prev] = NotNan::new(0.5).unwrap();
                                            pair_row[idx_current] = NotNan::new(0.0).unwrap();
                                            new_combinations.push(pair_row);
                                        }
                                    }
                                }
                            }

                            //track the current haplotype for pair generation
                            previous_haplotypes.push(ident_h.clone());
                        }
                    }
                }

                result.extend(new_combinations);
            }
        }
    }

    result
}

pub fn extend_resulting_table(
    representatives: &Vec<Haplotype>,
    event_posteriors: &Vec<(HaplotypeFractions, LogProb)>,
    prior: &PriorTypes,
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
            &prior,
            identical_haplotypes_map,
        );

        for combination in all_combinations {
            new_event_posteriors.push((HaplotypeFractions(combination), logprob.clone()));
        }
    }

    //convert logprobs to normal probabilities using
    let normal_probs: Vec<f64> = new_event_posteriors
        .iter()
        .map(|(_, logprob)| logprob.exp())
        .collect();

    let prob_sum: f64 = normal_probs.iter().sum();

    if prob_sum > 0.0 {
        //to avoid issues with dividing to 0
        //normalize probabilities and convert back to log space
        for (i, (_, logprob)) in new_event_posteriors.iter_mut().enumerate() {
            let normalized_prob = normal_probs[i] / prob_sum;
            *logprob = LogProb::from(Prob(normalized_prob));
        }
    }
    Ok((new_event_posteriors, all_haplotypes))
}
