use anyhow::Result;

use derive_builder::Builder;

use polars::frame::DataFrame;
use polars::prelude::*;

use crate::candidates::hla::{self};
use crate::candidates::hla::alignment;
use crate::candidates::hla::find_variants_from_cigar;
use crate::candidates::virus::sarscov2::write_to_vcf;

use bio::io::fasta;

use ndarray::Array2;
use petgraph::dot::{Config, Dot};
use petgraph::graph::{Graph, NodeIndex};
use petgraph::prelude::Dfs;

use rust_htslib::bcf::{header::Header, record::GenotypeAllele, Format, Writer};
use rust_htslib::faidx;
use seq_io::fasta::Record as OtherRecord;
use serde::Deserialize;

use std::collections::BTreeMap;
use std::collections::HashMap;

use std::fs;
use std::fs::File;
use std::io;

use std::io::Write;
use std::iter::FromIterator;
use std::path::PathBuf;

use std::process::Command;

#[allow(dead_code)]
#[derive(Builder, Clone)]
pub struct Caller {
    genome: PathBuf,
    lineages: PathBuf,
    output: PathBuf,
    threads: String,
}
impl Caller {
    pub fn call(&mut self) -> Result<()> {

        //align and sort
        alignment(
            &self.genome,
            &self.lineages,
            &self.threads,
            false,
            &self.output,
        )?;

        //find variants from cigar
        let (mut genotype_df, mut loci_df) = find_variants_from_cigar(
            &self.genome,
            &self.output.join("alignment_sorted.sam"),
        )
        .unwrap();


        //write locus-wise vcf files.
        write_to_vcf(&self.output, genotype_df, loci_df)?;

        Ok(())
    }
}