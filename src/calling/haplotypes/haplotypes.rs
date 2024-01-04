use crate::model::{AlleleFreq, Data, HaplotypeFractions, Likelihood, Marginal, Posterior, Prior};
use anyhow::Result;
use bio::stats::{bayesian::model::Model, probs::LogProb, PHREDProb, Prob};
use bv::BitVec;
use core::cmp::Ordering;
use derefable::Derefable;
use derive_builder::Builder;
use derive_deref::DerefMut;

use ordered_float::NotNan;
use ordered_float::OrderedFloat;

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
                && (&prob_absent_prob <= &Prob(0.05) || &prob_absent_prob >= &Prob(0.95))
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
pub(crate) struct HaplotypeVariants(
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
