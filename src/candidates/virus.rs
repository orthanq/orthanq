use anyhow::Result;

use bio_types::genome::AbstractInterval;
use csv::Reader as CsvReader;
use derive_builder::Builder;
use ndarray::Array2;
use ordered_float::NotNan;
use polars::{df, frame::DataFrame, prelude::NamedFrom, series::Series};
use rust_htslib::bcf::{header::Header, record::GenotypeAllele, Format, Writer};
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::record::Cigar, bam::Read, faidx};
use serde::Deserialize;
use std::collections::{BTreeMap, HashMap};
use std::convert::TryInto;

use std::iter::FromIterator;

use std::fs;
use std::path::PathBuf;
use std::process::Command;

use crate::candidates::hla;

#[allow(dead_code)]
#[derive(Builder, Clone)]
pub struct Caller {
    alleles: PathBuf,
    genome: PathBuf,
    output: Option<PathBuf>,
}
impl Caller {
    pub fn call(&self) -> Result<()> {
        //align and sort
        hla::alignment(&self.genome, &self.alleles)?;

        //find variants
        let (mut genotype_df, mut loci_df) =
            hla::find_variants_from_cigar(&self.genome, &"alignment_sorted.sam").unwrap();

        //write to vcf
        self.write_to_vcf(&genotype_df, &loci_df)?;

        Ok(())
    }

    fn write_to_vcf(&self, variant_table: &DataFrame, loci_table: &DataFrame) -> Result<()> {
        //Create VCF header
        let mut header = Header::new();
        //push contig names to the header.
        header.push_record(br#"##contig=<ID=1>"#);
        header.push_record(br#"##contig=<ID=6>"#);
        header.push_record(br#"##contig=<ID=7>"#);
        header.push_record(br#"##contig=<ID=8>"#);
        header.push_record(br#"##contig=<ID=9>"#);
        header.push_record(br#"##contig=<ID=11>"#);
        header.push_record(br#"##contig=<ID=16>"#);
        header.push_record(br#"##contig=<ID=X>"#);

        //push field names to the header.
        let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Variant is present in the haplotype (1) or not (0).">"#;
        header.push_record(header_gt_line.as_bytes());
        let header_locus_line = r#"##FORMAT=<ID=C,Number=1,Type=Integer,Description="Locus is covered by the haplotype (1) or not (0).">"#;
        header.push_record(header_locus_line.as_bytes());

        //push sample names into the header
        for sample_name in variant_table.get_column_names().iter().skip(2) {
            header.push_sample(sample_name.as_bytes());
        }
        fs::create_dir_all(self.output.as_ref().unwrap())?;
        let mut vcf = Writer::from_path(
            format!(
                "{}.vcf",
                self.output
                    .as_ref()
                    .unwrap()
                    .join("candidates.vcf")
                    .display()
            ),
            &header,
            true,
            Format::Vcf,
        )
        .unwrap();

        let _id_iter = variant_table["ID"].i64().unwrap().into_iter();
        for row_index in 0..variant_table.height() {
            let mut record = vcf.empty_record();
            let mut variant_iter = variant_table["Index"].utf8().unwrap().into_iter();
            let mut id_iter = variant_table["ID"].i64().unwrap().into_iter();
            let splitted = variant_iter
                .nth(row_index)
                .unwrap()
                .unwrap()
                .split(',')
                .collect::<Vec<&str>>();
            let chrom = splitted[0];
            let pos = splitted[1];
            let id = id_iter.nth(row_index).unwrap().unwrap();
            let ref_base = splitted[2];
            let alt_base = splitted[3];
            let alleles: &[&[u8]] = &[ref_base.as_bytes(), alt_base.as_bytes()];
            let rid = vcf.header().name2rid(chrom.as_bytes()).unwrap();

            record.set_rid(Some(rid));
            record.set_pos(pos.parse::<i64>().unwrap() - 1);
            record.set_id(id.to_string().as_bytes()).unwrap();
            record.set_alleles(alleles).expect("Failed to set alleles");

            //push genotypes
            let mut all_gt = Vec::new();
            for column_index in 2..variant_table.width() {
                let gt = variant_table[column_index]
                    .i32()
                    .unwrap()
                    .into_iter()
                    .nth(row_index)
                    .unwrap()
                    .unwrap();
                let gt = GenotypeAllele::Phased(gt);
                //vg requires the candidate variants phased, so make 0 -> 0|0 and 1 -> 1|1
                //side note about how push_genotypes() works:
                //we push the genotypes two times into the one dimensional vector because
                //push_genotypes() expects a 'flattened' two dimensional vector
                //thereby pushing into the vector two times also means the same thing for push_genotypes() in the end
                all_gt.push(gt);
                all_gt.push(gt);
            }
            record.push_genotypes(&all_gt).unwrap();

            //push loci
            let mut all_c = Vec::new();
            for column_index in 2..loci_table.width() {
                //it doesnt have the ID column so that it starts from 1
                let c = loci_table[column_index]
                    .i32()
                    .unwrap()
                    .into_iter()
                    .nth(row_index)
                    .unwrap()
                    .unwrap();
                all_c.push(c);
            }
            record.push_format_integer(b"C", &all_c)?;
            vcf.write(&record).unwrap();
        }
        Ok(())
    }
}
