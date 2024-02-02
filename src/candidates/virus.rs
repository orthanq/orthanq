use anyhow::Result;

use derive_builder::Builder;

use polars::frame::DataFrame;
use polars::prelude::*;
use polars::df;
use rust_htslib::bcf::{header::Header, record::GenotypeAllele, Format, Writer};
// use rust_htslib::{bam::Read};
use std::ffi::OsStr;
use std::fs;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::process;
use std::io::Write;
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
        //prepare genome and alleles for the virus (currently, only sars-cov2) (ncbi-datasets-cli=16.4.4 package should be in the requiremenbts)
        //first, create the output dir for database setup
        fs::create_dir_all(self.output.as_ref().unwrap())?;

        //output dir for metadata
        let virus_metadata_dir = self.output.as_ref().unwrap().join("metadata.tsv");

        //download and prepare required metadata as tsv
        //datasets summary virus genome taxon SARS-CoV-2 --as-json-lines | dataformat tsv virus-genome --fields accession,completeness,release-date,virus-pangolin > sarscov2.tsv
        
        //first, download metada
        let download_metadata = {
            Command::new("datasets")
                .arg("summary")
                .arg("virus")
                .arg("genome")
                .arg("taxon")
                .arg("SARS-CoV-2")
                .arg("--as-json-lines")
                .stdout(Stdio::piped())
                .spawn()
                .expect("failed to execute the vg giraffe process")
        };
        println!("virus metada download is complete.");

        //then, pipe the downloaded jsonl to tsv with dataformat
        let process_into_tsv = Command::new("dataformat")
            .arg("tsv")
            .arg("virus-genome")
            .stdin(download_metadata.stdout.unwrap())
            .arg("--fields")
            .arg("accession,completeness,release-date,virus-pangolin")
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();

        //write the prepared tsv to file
        let output = process_into_tsv
            .wait_with_output()
            .expect("failed to wait on child");
        let mut f = std::fs::File::create(virus_metadata_dir.clone())?;
        f.write_all(&output.stdout)?;
        println!("virus metada preparation is complete.");

        //comment out the next line
        // let virus_metadata_dir = &"metadata.tsv";

        //retrieve accession ids from each lineage with the oldest submission and complete genomes.
        let virus_metadata_df = CsvReader::from_path(virus_metadata_dir)?.with_delimiter(b'\t').has_header(true).finish()?;
        dbg!(&virus_metadata_df);
        
        //filter for complete genomes and non-null pangolin classification

        // create a mask to filter out null values
        let mask = virus_metadata_df.column("Virus Pangolin Classification")?.is_not_null();

        // apply the filter on a DataFrame
        let filtered_df_annotated = virus_metadata_df.filter(&mask)?;
        dbg!(&filtered_df_annotated);

        //then, only include complete genomes
        let filter_complete = filtered_df_annotated["Completeness"].utf8().unwrap().into_iter().map(|x| 
            if x.unwrap() == "COMPLETE" { true } else { false }).collect();
        let filtered_df_complete = filtered_df_annotated.filter(&filter_complete).unwrap();
        dbg!(&filtered_df_complete);

        //group by pangolin lineage and sort by release date

        // ordering of the columns
        let descending = vec![false, false];
        // columns to sort by
        let by = &["Virus Pangolin Classification", "Release date"];
        // do the sort operation
        let sorted = filtered_df_complete.sort(by, descending)?;
        dbg!(&sorted);

        //group by pangolin lineage and get the first entry (which is the oldest)
        let first = sorted.groupby(["Virus Pangolin Classification"])?.select(["Accession"]).first()?;
        dbg!(&first);

        //download genomes by the accession
        //e.g. datasets download virus genome accession OY726946.1,OY299705.1
        dbg!(&first["Accession_first"]);
        let accessions_opt: Vec<Option<&str>> = first["Accession_first"].utf8()?.into_iter().collect();
        let accessions: Vec<&str> = accessions_opt.iter().map(|a|a.unwrap()).collect();
        dbg!(&accessions);
        dbg!(&accessions.len());

        //write accessions to file
        let accessions_input = self.output.as_ref().unwrap().join("accessions.csv");
        let mut wtr = csv::Writer::from_path(&accessions_input)?;
        for accession in accessions.iter(){
            wtr.write_record(vec![accession.clone()])?;
        }
        wtr.flush(); //make sure everything is sucsessfully written to the file 

        //download genomes given by accession
        let download_genomes = {
            Command::new("datasets")
                .arg("download")
                .arg("virus")
                .arg("genome")
                .arg("accession")
                .arg("--inputfile")
                .arg(accessions_input)
                .arg("--debug")
                .status()
                .expect("failed to download sequences")
        };
        println!("The download of sars-cov-2 genomes with given accession ids was exited with: {}", download_genomes);

        process::exit(0x0100);

        //align and sort
        hla::alignment(&self.genome, &self.alleles)?;

        //find variants
        let (genotype_df, loci_df) =
            hla::find_variants_from_cigar(&self.genome, &"alignment_sorted.sam").unwrap();

        //write to vcf
        self.write_to_vcf(&genotype_df, &loci_df)?;

        Ok(())
    }

    fn write_to_vcf(&self, variant_table: &DataFrame, loci_table: &DataFrame) -> Result<()> {
        // dbg!(&variant_table);
        //Create VCF header
        let mut header = Header::new();
        //push contig names to the header.
        header.push_record(br#"##contig=<ID=NC_045512.2>"#);

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
        let mut vcf = Writer::from_path(
            format!(
                "{}.vcf",
                self.output.as_ref().unwrap().join("candidates").display()
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
