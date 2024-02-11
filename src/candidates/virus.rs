use anyhow::Result;

use derive_builder::Builder;

use polars::frame::DataFrame;
use polars::prelude::*;

use rust_htslib::bcf::{header::Header, record::GenotypeAllele, Format, Writer};
use std::fs;

use std::path::PathBuf;
use std::process::{Command, Stdio};

use crate::candidates::hla;

use std::io::Write;

#[allow(dead_code)]
#[derive(Builder, Clone)]
pub struct Caller {
    output: Option<PathBuf>,
    threads: String,
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
        println!("virus metadata preparation is complete.");

        //comment out the next line for testing purposees (above process can take ~6 hours)
        // let virus_metadata_dir = "/home/uzuner/Documents/backup_orthanq_metadata/metadata.tsv";

        //retrieve accession ids from each lineage with the oldest submission and complete genomes.
        let virus_metadata_df = CsvReader::from_path(virus_metadata_dir)?
            .with_delimiter(b'\t')
            .has_header(true)
            .finish()?;
        dbg!(&virus_metadata_df);

        //filter for complete genomes and non-null pangolin classification

        //then, only include complete genomes
        let filter_complete = virus_metadata_df["Completeness"]
            .utf8()
            .unwrap()
            .into_iter()
            .map(|x| {
                if x.unwrap() == "COMPLETE" {
                    true
                } else {
                    false
                }
            })
            .collect();
        let filtered_df_complete = virus_metadata_df.filter(&filter_complete).unwrap();
        dbg!(&filtered_df_complete);

        // create a mask to filter out null values
        let mask = filtered_df_complete
            .column("Virus Pangolin Classification")?
            .is_not_null();

        // apply the filter on a DataFrame
        let filtered_df_annotated = filtered_df_complete.filter(&mask)?;
        dbg!(&filtered_df_annotated);

        //group by pangolin lineage and sort by release date

        // ordering of the columns
        let descending = vec![false, false];
        // columns to sort by
        let by = &["Virus Pangolin Classification", "Release date"];
        // do the sort operation
        let sorted = filtered_df_annotated.sort(by, descending)?;

        dbg!(&sorted);

        //group by pangolin lineage and get the first entry (which is the oldest)
        let grouped_df = sorted
            .groupby(["Virus Pangolin Classification"])?
            .select(["Accession"])
            .first()?;
        dbg!(&grouped_df);

        //download genomes by the accession
        //e.g. datasets download virus genome accession OY726946.1,OY299705.1
        dbg!(&grouped_df["Accession_first"]);
        let accessions_opt: Vec<Option<&str>> =
            grouped_df["Accession_first"].utf8()?.into_iter().collect();
        let accessions: Vec<&str> = accessions_opt.iter().map(|a| a.unwrap()).collect();
        // dbg!(&accessions);
        dbg!(&accessions.len());

        //write accessions to file
        let accessions_input = self.output.as_ref().unwrap().join("accessions.csv");
        let mut wtr = csv::Writer::from_path(&accessions_input)?;
        for accession in accessions.iter() {
            wtr.write_record(vec![accession])?;
        }
        wtr.flush()?; //make sure everything is sucsessfully written to the file

        //download genomes given by accession, this will create a file called ncbi_datasets.zip
        let data_package_lineages = self
            .output
            .as_ref()
            .unwrap()
            .join("ncbi_datasets_lineages.zip");

        let download_genomes = {
            Command::new("datasets")
                .arg("download")
                .arg("virus")
                .arg("genome")
                .arg("accession")
                .arg("--inputfile")
                .arg(accessions_input)
                .arg("--filename")
                .arg(&data_package_lineages)
                .arg("--debug")
                .status()
                .expect("failed to download sequences")
        };
        println!(
            "The download of sars-cov-2 genomes with given accession ids was exited with: {}",
            download_genomes
        );

        //unzip the downloaded package
        let unzip_package = {
            Command::new("unzip")
                .arg(data_package_lineages)
                .arg("-d")
                .arg(self.output.as_ref().unwrap().join("lineages"))
                .status()
                .expect("failed to unzip ncbi data package")
        };
        println!("All genomes data package was unzipped: {}", unzip_package);

        //download the oldest submitted sars-cov-2 genome to use as reference
        // first, sort by release
        let sorted = filtered_df_complete.sort(&["Release date"], false)?;
        dbg!(&sorted);

        //write to csv
        // let mut output_file: File = File::create("out.csv").unwrap();
        // CsvWriter::new(&mut output_file)
        //     .has_header(true)
        //     .finish(&mut sorted)
        //     .unwrap();

        //download the oldest submitted genome to be used as reference genome
        // let oldest_accession = sorted.select(["Accession"]).iter();
        let accessions = sorted.column("Accession")?.utf8()?;
        let accessions_to_vec: Vec<Option<&str>> = accessions.into_iter().collect();
        let oldest_accession = accessions_to_vec[0];

        //download the genome
        let data_package_reference = self
            .output
            .as_ref()
            .unwrap()
            .join("ncbi_datasets_reference_genome.zip");

        let download_reference_genome = {
            Command::new("datasets")
                .arg("download")
                .arg("virus")
                .arg("genome")
                .arg("accession")
                .arg(oldest_accession.unwrap())
                .arg("--filename")
                .arg(&data_package_reference)
                .arg("--debug")
                .status()
                .expect("failed to download sequences")
        };
        println!(
            "The download of reference genome for sars-cov-2 was exited with: {}",
            download_reference_genome
        );

        //unzip the downloaded package
        let unzip_package = {
            Command::new("unzip")
                .arg(data_package_reference)
                .arg("-d")
                .arg(self.output.as_ref().unwrap().join("reference_genome"))
                .status()
                .expect("failed to unzip ncbi data package for reference genome")
        };
        println!(
            "The reference genome data package was unzipped: {}",
            unzip_package
        );

        //then, here are the dirs for reference genome and all lineages
        let lineages_dir = self
            .output
            .as_ref()
            .unwrap()
            .join("lineages/ncbi_dataset/data/genomic.fna");
        let reference_genome_dir = self
            .output
            .as_ref()
            .unwrap()
            .join("reference_genome/ncbi_dataset/data/genomic.fna");

        //index reference genome
        let faidx = {
            Command::new("samtools")
                .arg("faidx")
                .arg(&reference_genome_dir)
                .status()
                .expect("failed to unzip ncbi data package for reference genome")
        };
        println!(
            "Samtools faidx for reference genome was exited with: {}",
            faidx
        );

        //align and sort

        hla::alignment(
            &reference_genome_dir,
            &lineages_dir,
            &self.threads,
            false,
            self.output.as_ref().unwrap(),
        )?;

        //find variants
        let (genotype_df, loci_df) = hla::find_variants_from_cigar(
            &reference_genome_dir,
            &self.output.as_ref().unwrap().join("alignment_sorted.sam"),
        )
        .unwrap();

        //replace accession ids with the corresponding lineage names
        let genotype_lineages = accession_to_lineage(&genotype_df, &grouped_df)?;
        let loci_lineages = accession_to_lineage(&loci_df, &grouped_df)?;

        //write to vcf
        self.write_to_vcf(&genotype_lineages, &loci_lineages)?;

        Ok(())
    }

    fn write_to_vcf(&self, variant_table: &DataFrame, loci_table: &DataFrame) -> Result<()> {
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

//accession_to_lineage simply converts accession ids that represent each lineage to lineage names as well as nextstrain clade
fn accession_to_lineage(df: &DataFrame, accession_to_lineage_df: &DataFrame) -> Result<DataFrame> {
    dbg!(&df);
    dbg!(&accession_to_lineage_df);

    //read lineage to clade resource
    let cargo_dir = env!("CARGO_MANIFEST_DIR");
    let file_dir = format!("{}/resources/clade_to_lineages/cladeToLineages.tsv", cargo_dir);
    let lin_to_clade = CsvReader::from_path(&file_dir)?
        .with_delimiter(b'\t')
        .has_header(true)
        .finish()?;
    dbg!(&lin_to_clade);

    //clone the input df
    let mut renamed_df = df.clone();

    //prepare iterables of columns beforehand
    let mut iter_pangolin = accession_to_lineage_df["Virus Pangolin Classification"]
        .utf8()
        .unwrap()
        .into_iter();
    let mut iter_accession = accession_to_lineage_df["Accession_first"]
        .utf8()
        .unwrap()
        .into_iter();
    let mut lineages_in_resource = lin_to_clade["lineage"].utf8().unwrap();
    let mut clades_in_resource = lin_to_clade["clade"].utf8().unwrap();

    //rename accessions to lineages and if present, clades
    iter_pangolin
        .into_iter()
        .zip(iter_accession)
        .for_each(|(png, acc)| {
            let mut rename_str = format!("no clade ({})", png.unwrap());
            lineages_in_resource
                .into_iter()
                .zip(clades_in_resource.into_iter())
                .for_each(|(lng, cld)| {
                    if png == lng {
                        rename_str = format!("{} ({})", cld.unwrap(), lng.unwrap());
                    }
                });
            renamed_df.rename(acc.unwrap(), &rename_str);
        });
    dbg!(&renamed_df);
    Ok(renamed_df)
}
