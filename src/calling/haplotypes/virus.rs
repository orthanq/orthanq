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
    outcsv: Option<PathBuf>,
    prior: String,
}

impl Caller {
    pub fn call(&mut self) -> Result<()> {
        Ok(())
    }
}
