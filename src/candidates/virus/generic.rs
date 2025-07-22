use anyhow::Result;

use derive_builder::Builder;

use crate::candidates::hla::alignment;
use crate::candidates::hla::find_variants_from_cigar;

use std::io::Write;

use std::path::PathBuf;

use polars::frame::DataFrame;
use polars::prelude::*;

use rust_htslib::bcf::{header::Header, record::GenotypeAllele, Format, Writer};

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
            &"virus",
            &self.genome,
            &self.lineages,
            &self.threads,
            false,
            &self.output,
        )?;

        //find variants from cigar
        let (genotype_df, loci_df) =
            find_variants_from_cigar(&self.genome, &self.output.join("viruses_alignment_sorted.bam"))
                .unwrap();

        //write locus-wise vcf files.
        write_to_vcf(&self.output, genotype_df, loci_df)?;

        Ok(())
    }
}

pub fn write_to_vcf(
    outdir: &PathBuf,
    variant_table: DataFrame,
    loci_table: DataFrame,
) -> Result<()> {
    // dbg!(&variant_table);
    //Create VCF header
    let mut header = Header::new();

    //get contig name of the reference
    let first_row_index = variant_table["Index"]
        .utf8()
        .unwrap()
        .into_iter()
        .nth(0)
        .unwrap()
        .unwrap()
        .split(',')
        .collect::<Vec<&str>>();

    //push contig name to the header
    let header_contig_line = format!(r#"##contig=<ID={}>"#, first_row_index[0]);
    header.push_record(header_contig_line.as_bytes());

    //push field names to the header.
    let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Variant is present in the haplotype (1) or not (0).">"#;
    header.push_record(header_gt_line.as_bytes());
    let header_locus_line = r#"##FORMAT=<ID=C,Number=1,Type=Integer,Description="Locus is covered by the haplotype (1) or not (0).">"#;
    header.push_record(header_locus_line.as_bytes());

    //push sample names into the header
    for sample_name in variant_table.get_column_names().iter().skip(2) {
        header.push_sample(sample_name.as_bytes());
    }
    // fs::create_dir_all(self.output.as_ref().unwrap())?;
    // let outdir: &_ = outdir;

    let mut vcf = Writer::from_path(
        format!("{}.vcf", outdir.join("candidates").display()),
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
