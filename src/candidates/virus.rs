use anyhow::Result;

use derive_builder::Builder;

use polars::frame::DataFrame;
use polars::prelude::*;

use crate::candidates::hla;
use bio::io::fasta::{Record, Writer as FastaWriter};
use reqwest;
use rust_htslib::bcf::{header::Header, record::GenotypeAllele, Format, Writer};
use rust_htslib::faidx;
use seq_io::fasta::{Reader, Record as OtherRecord};
use std::collections::HashSet;
use std::convert::TryInto;
use std::error::Error;
use std::fs;
use std::io;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::process::ExitStatus;
use std::process::{Command, Stdio};

#[allow(dead_code)]
#[derive(Builder, Clone)]
pub struct Caller {
    output: Option<PathBuf>,
    threads: String,
}
impl Caller {
    pub fn call(&self) -> Result<()> {
        //prepare genome and alleles for the virus (currently, only sars-cov2)

        //define cargo dir and state the location of cladetolineages file.
        let cargo_dir = env!("CARGO_MANIFEST_DIR");
        let file_path_clade_to_lineages = format!(
            "{}/resources/clade_to_lineages/cladeToLineages.tsv",
            cargo_dir
        );
        dbg!(&file_path_clade_to_lineages);

        //first, create the output dir for database setup
        fs::create_dir_all(self.output.as_ref().unwrap())?;

        //access sarscov2 data package via REST API of NCBI
        let path_to_ncbi_dataset = self.output.as_ref().unwrap().join("ncbi_dataset.zip");

        //the sarscov2 genome data package is downloaded and written to file.

        // // --START-- ALTERNATIVE way to download data package (very very slow, might take days): download accessions separately and combine ( bc. downloading the whole sarscov2 data package at once might run into timeout errors on the server side.)

        // // //output dir for metadata
        // // let virus_metadata_path = self.output.as_ref().unwrap().join("metadata.tsv");

        // // let download_metadata = {
        // //     Command::new("datasets")
        // //         .arg("summary")
        // //         .arg("virus")
        // //         .arg("genome")
        // //         .arg("taxon")
        // //         .arg("SARS-CoV-2")
        // //         .arg("--as-json-lines")
        // //         .stdout(Stdio::piped())
        // //         .spawn()
        // //         .expect("failed to execute the vg giraffe process")
        // // };
        // // println!("virus metada download is complete.");

        // // //then, pipe the downloaded jsonl to tsv with dataformat
        // // let process_into_tsv = Command::new("dataformat")
        // //     .arg("tsv")
        // //     .arg("virus-genome")
        // //     .stdin(download_metadata.stdout.unwrap())
        // //     .arg("--fields")
        // //     .arg("accession")
        // //     .stdout(Stdio::piped())
        // //     .spawn()
        // //     .unwrap();

        // // //write the prepared tsv to file
        // // let output = process_into_tsv
        // //     .wait_with_output()
        // //     .expect("failed to wait on child");
        // // let mut f = std::fs::File::create(virus_metadata_path.clone())?;
        // // f.write_all(&output.stdout)?;
        // // println!("virus metadata preparation is complete.");

        // // //comment out the next line for testing purposees (above process can take ~6 hours)
        // // let virus_metadata_path = "/home/hamdiyeuzuner/Documents/test_orthanq/metadata.tsv";

        // // //grab accession ids and download them one by one

        // // let accessions_df = CsvReader::from_path(virus_metadata_path)?
        // //     .with_delimiter(b'\t')
        // //     .has_header(true)
        // //     .finish()?;
        // // dbg!(&accessions_df);
        // // let accessions_vec: Vec<&str> = accessions_df[0]
        // //     .utf8()
        // //     .unwrap()
        // //     .into_iter()
        // //     .map(|a|a.unwrap())
        // //     .collect();
        // // dbg!(&accessions_df);

        // // let len_accessions = accessions_vec.len();
        // // dbg!(&len_accessions);
        // // println!("genomes to download: {}", len_accessions.to_string());

        // // //monitor the progress of download with the counter
        // // let mut counter_progress = 0;

        // // for accession in accessions_vec.iter(){
        // //     dbg!(&accession);
        // //     //create the folder to collect all genomes
        // //     let genomes_dir = self.output.as_ref().unwrap().join("genomes");
        // //     fs::create_dir_all(genomes_dir)?;

        // //     let path_to_genome = self.output.as_ref().unwrap().join(format!("genomes/{}.zip", accession.to_string()));

        // //     let download_genomes_status = download_accession(accession.to_string(), &path_to_genome).unwrap();
        // //     // println!("status: {}", &download_genomes_status); //show status

        // //     //if the download has success, print the progress info, if not try once again and if the process is still not successful error, then panic
        // //     if download_genomes_status.success() {
        // //         counter_progress += 1;
        // //         println!("{}/{} genome downloaded", counter_progress.to_string(),len_accessions.to_string())
        // //     } else {
        // //         println!("genome could not be downloaded, retrying");
        // //         let download_genomes_status = download_accession(accession.to_string(), &path_to_genome).unwrap();
        // //         println!("status: {}", &download_genomes_status); //show status
        // //         if !download_genomes_status.success() {
        // //             panic!("genome could not be downloaded");
        // //         }
        // //     }
        // // }

        // // println!("virus genome data download is complete.");

        // // //download each genome by accession
        // // let download_genomes = {
        // //     Command::new("datasets")
        // //         .arg("download")
        // //         .arg("virus")
        // //         .arg("genome")
        // //         .arg("accession")
        // //         .arg("--inputfile")
        // //         .arg(virus_metadata_path)
        // //         .status()
        // //         .expect("failed to execute the download process")
        // // };
        // // println!("virus genome data download is complete.");
        // //todo: combine all genomes to have the sars-cov-2 data package in the end.
        // // --END-- ALTERNATIVE way

        /////////////////////////////////////////////////////////
        //1-) Download the whole genome package for SARSCoV-2. //
        /////////////////////////////////////////////////////////

        //Download the whole genome package at once (occasionally might run into timeout errors) (much faster than the alternative approach, not more half an hour usually)
        if !path_to_ncbi_dataset.exists() {
            let mut child = {
                Command::new("datasets")
                    .arg("download")
                    .arg("virus")
                    .arg("genome")
                    .arg("taxon")
                    .arg("SARS-CoV-2")
                    .arg("--complete-only")
                    .arg("--filename")
                    .arg(&path_to_ncbi_dataset)
                    .status()
                    .expect("failed to download sequences")
            };
            println!("SARS-CoV-2 data package was downloaded.")
        }

        //unzip the downloaded package
        let unzip_package = {
            Command::new("unzip")
                .arg(path_to_ncbi_dataset)
                .arg("-d")
                .arg(self.output.as_ref().unwrap().join("ncbi_genomes_all"))
                .status()
                .expect("failed to unzip ncbi data package")
        };
        println!("All genomes data package was unzipped: {}", unzip_package);

        /////////////////////////////////////////////////////////////////////////////////////////
        //2-) Perform quality control on the sequences & write surviving accession IDs to file.//
        /////////////////////////////////////////////////////////////////////////////////////////

        //select only fasta records containing less than a certain threshold for ambiguous bases, threshold = 0.05
        //for that, count the number of ambiguous bases i.e. N, R, etc. and divide them by the length of the fasta
        //also, for a check of genome completeness, check LR877184.1, only include genome lenghts with greater than 29k ref:  (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7161481/)

        //location of fasta belonging to all genomes.
        let all_genomes = self
            .output
            .as_ref()
            .unwrap()
            .join("ncbi_genomes_all/ncbi_dataset/data/genomic.fna");

        //initialize the vector to collect quality checked accession ids
        let mut quality_accessions: Vec<String> = Vec::new();

        //find the records that are above the threshold and collect accessions to vec.
        //iterative operation over the reference genome is done with seq_io crate, instead of rust-htslib, because iterating over
        //reference genome that is ~80GB leads to a memory overflow and system crash.

        let mut reader = seq_io::fasta::Reader::from_path(&all_genomes).unwrap();
        println!("below process takes around half an hour.");
        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");

            //find the name, length and sequence of the fasta record
            let l_name = record.id()?;
            let l_len = record.seq_lines().fold(0, |l, seq| l + seq.len());
            let l_seq = String::from_utf8(record.seq().to_vec()).unwrap();

            //find the number of ambiguous bases in the record
            let n_ambiguous = l_seq
                .chars()
                .filter(|x| !vec!["A", "G", "T", "C"].contains(&x.to_string().as_str()))
                .count();

            //calculate the ratio of ambiguous bases in the fasta record
            let ratio_of_ambiguous = n_ambiguous as f32 / l_len as f32;

            //write the record to file if the ratio exceeds the defined threshold
            let threshold_ambiguous: f32 = 0.05; //todo: must be configured later

            if (ratio_of_ambiguous < threshold_ambiguous) && (l_len > 29000) {
                //second check is necessary for double checking genome completeness
                quality_accessions.push(l_name.to_string());
            }
        }
        println!("above process done.");

        //write accessions to file
        let accessions_input = self.output.as_ref().unwrap().join("accessions.csv");
        let mut wtr = csv::Writer::from_path(&accessions_input)?;
        for accession in quality_accessions.iter() {
            wtr.write_record(vec![accession])?;
        }
        wtr.flush()?; //make sure everything is sucsessfully written to the file

        // //uncomment the next line, it's only for testing
        // let accessions_input = self.output.as_ref().unwrap().join("accessions.csv");

        ///////////////////////////////////////////////////////////////////////////////////////
        //3-) Download metadata only for the accessions coming from the previous step.         //
        ///////////////////////////////////////////////////////////////////////////////////////

        println!("virus metada download started.");
        let download_metadata = {
            Command::new("datasets")
                .arg("summary")
                .arg("virus")
                .arg("genome")
                .arg("accession")
                .arg("--inputfile")
                .arg(accessions_input)
                .arg("--as-json-lines")
                .stdout(Stdio::piped())
                .spawn()
                .expect("failed to execute the vg giraffe process")
        };
        println!("virus metada download is complete.");

        //then, pipe the downloaded jsonl to tsv with dataformat
        println!("virus metadata preparation started.");

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

        //output dir for metadata
        let virus_metadata_path = self.output.as_ref().unwrap().join("metadata.tsv");

        let output = process_into_tsv
            .wait_with_output()
            .expect("failed to wait on child");
        let mut f = std::fs::File::create(virus_metadata_path.clone())?;
        f.write_all(&output.stdout)?;
        println!("virus metadata preparation is complete.");

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //4-) Find representative accession IDs for each Nextstrain clade by taking the oldest submitted genome.    //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // //retrieve accession ids from each clade with the oldest submission.
        let virus_metadata_df = CsvReader::from_path(virus_metadata_path)?
            .with_delimiter(b'\t')
            .has_header(true)
            .finish()?;
        dbg!(&virus_metadata_df.shape());

        //merge with cladeToLineages.tsv to group for nextstrain clades
        //read lineage to clade resource
        let lin_to_clade = CsvReader::from_path(&file_path_clade_to_lineages)?
            .with_delimiter(b'\t')
            .has_header(true)
            .finish()?;
        dbg!(&lin_to_clade.shape());

        //merge virus_metadata_df with lin_to_clade
        let virus_metadata_wnextstrain = virus_metadata_df.inner_join(
            &lin_to_clade,
            ["Virus Pangolin Classification"],
            ["lineage"],
        )?;
        dbg!(&virus_metadata_wnextstrain.shape());

        //then, only include complete genomes
        let filter_complete = virus_metadata_wnextstrain["Completeness"]
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
        let filtered_df_complete = virus_metadata_wnextstrain.filter(&filter_complete).unwrap();
        dbg!(&filtered_df_complete.shape());

        //group by clade and sort by release date

        // ordering of the columns
        let descending = vec![false, false];
        // columns to sort by
        let by = &["clade", "Release date"];
        // do the sort operation
        let mut sorted = filtered_df_complete.sort(by, descending)?;
        dbg!(&sorted.shape());

        //for debuggin, write merged_sorted_df to file
        let merged_sorted_path = self.output.as_ref().unwrap().join("merged_sorted.csv");

        let mut file = std::fs::File::create(merged_sorted_path).expect("could not create file");

        CsvWriter::new(&mut file)
            .has_header(true)
            .with_delimiter(b',')
            .finish(&mut sorted);

        //group by clade and get the first entry (which is the oldest)
        let grouped_df = sorted.groupby(["clade"])?.select(["Accession"]).first()?;
        dbg!(&grouped_df.shape());

        //collect first accessions to a vector
        let selected_accessions_opt: Vec<Option<&str>> =
            grouped_df["Accession_first"].utf8()?.into_iter().collect();

        //collect the accessions in a HashSet to avoid the same accessions that represent multiple clades to be added.
        //todo: think about how to include other lineages
        let selected_accessions: HashSet<_> =
            selected_accessions_opt.iter().map(|a| a.unwrap()).collect();
        dbg!(&selected_accessions.len());

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //5-) Write the representative accessions to csv and then take subset using seqtk and write to fasta (fast) //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        let selected_accessions_csv = self
            .output
            .as_ref()
            .unwrap()
            .join("selected_accessions.csv");
        let mut selected_accessions_wtr = csv::Writer::from_path(&selected_accessions_csv)?;
        for accession in selected_accessions.iter() {
            selected_accessions_wtr.write_record(vec![accession])?;
        }
        selected_accessions_wtr.flush()?; //make sure everything is sucsessfully written to the file

        //then take subset of the main fasta with seqtk (fast)
        let selected_accessions_fna = self
            .output
            .as_ref()
            .unwrap()
            .join("selected_quality_accessions.fna");

        //location of fasta belonging to all genomes.
        let all_genomes = self
            .output
            .as_ref()
            .unwrap()
            .join("ncbi_genomes_all/ncbi_dataset/data/genomic.fna");

        //take subset of selected accessions and write to file.
        println!("subsetting of fasta with quality lineages started");
        let seqtk_subset = {
            Command::new("seqtk")
                .arg("subseq")
                .arg(&all_genomes)
                .arg(selected_accessions_csv)
                .output()
                .expect("failed to take subset of FASTA")
        };
        println!("subsetting of fasta done.");

        let stdout = String::from_utf8(seqtk_subset.stdout).unwrap();
        fs::write(&selected_accessions_fna, stdout)
            .expect("Unable to write subsetted fasta to file");
        println!("subsetted fasta written to file. ");

        //////////////////////////////////////////////////////
        //6-) Find the reference genome and write to fasta. //
        //////////////////////////////////////////////////////

        let reference_fna = self.output.as_ref().unwrap().join("reference_genome.fna");

        //find the oldest submitted sars-cov-2 genome
        let sorted_df_by_release = filtered_df_complete.sort(&["Release date"], false)?;
        let oldest_genome_accession: String = sorted_df_by_release["Accession"]
            .utf8()?
            .into_iter()
            .collect::<Vec<_>>()
            .iter()
            .map(|a| a.unwrap().to_string())
            .collect::<Vec<_>>()[0]
            .to_string();

        //write to csv (file input is required for seqtk)
        let reference_csv = self.output.as_ref().unwrap().join("reference.csv");
        let mut reference_wtr = csv::Writer::from_path(&reference_csv)?;
        reference_wtr.write_record(vec![oldest_genome_accession.clone()])?;
        reference_wtr.flush()?; //make sure everything is sucsessfully written to the file

        //take subset of the reference genome from selected quality lineages fasta.
        println!("subsetting of fasta for the reference genome started");
        let seqtk_subset = {
            Command::new("seqtk")
                .arg("subseq")
                .arg(&selected_accessions_fna)
                .arg(reference_csv)
                .output()
                .expect("failed to take subset of FASTA")
        };
        println!("subsetting of fasta done.");

        let stdout = String::from_utf8(seqtk_subset.stdout).unwrap();
        fs::write(&reference_fna, stdout).expect("Unable to write subsetted fasta to file");
        println!("subsetted fasta written to file. ");

        //index reference genome
        let faidx = {
            Command::new("samtools")
                .arg("faidx")
                .arg(&reference_fna)
                .status()
                .expect("failed to unzip ncbi data package for reference genome")
        };
        println!(
            "Samtools faidx for reference genome was exited with: {}",
            faidx
        );

        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //7-) Aligned selected quality accessions to the reference genome and sort the resulting alignment file.//
        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        //DELETE AFTER TESTING ABOVE PART.
        // //following two lines must be commented out, they're only for testing.
        // let quality_lineages_path = self.output.as_ref().unwrap().join("quality_accessions.fna");
        // let reference_genome_path = self.output.as_ref().unwrap().join("reference_genome.fna");

        hla::alignment(
            &reference_fna,
            &selected_accessions_fna,
            &self.threads,
            false,
            self.output.as_ref().unwrap(),
        )?;

        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //8-) Find variants from the alignment, convert accessions to lineage names in the resulting dataframe. //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        let (genotype_df, loci_df) = hla::find_variants_from_cigar(
            &reference_fna,
            &self.output.as_ref().unwrap().join("alignment_sorted.sam"),
        )
        .unwrap();

        //replace accession ids with the corresponding lineage names
        //the performance of the following function usage can be improved for performance reasons (the whole virus_metadata_df might not have to be iterated over.)
        let genotype_lineages = accession_to_lineage(
            &genotype_df,
            &virus_metadata_df,
            &file_path_clade_to_lineages,
        )?;
        let loci_lineages =
            accession_to_lineage(&loci_df, &virus_metadata_df, &file_path_clade_to_lineages)?;

        ///////////////////////////////////////////////
        //9-) Write the candidate variants to fasta. //
        ///////////////////////////////////////////////

        self.write_to_vcf(&genotype_lineages, &loci_lineages)?;

        Ok(())
    }

    #[tokio::main]
    async fn get_data_ncbi(path_to_ncbi_dataset: &PathBuf) -> Result<()> {
        //try the REST API of NCBI
        // chaining .await will yield our query result
        let url = format!(
            "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/virus/taxon/2697049/genome/download?refseq_only=false&annotated_only=false&complete_only=true&include_sequence=GENOME&pangolin_classification=XBB.1.5.109",
        );
        //first access this: NC_045512.2
        // let url = format!(
        //     "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/virus/accession/NC_045512.2/genome/download?refseq_only=true&annotated_only=true&host=human&complete_only=true&include_sequence=GENOME&include_sequence=CDS&include_sequence=PROTEIN",
        // );
        let client = reqwest::Client::new();
        let response = client
            .get(url)
            // confirm the request using send()
            .send()
            .await
            // the rest is the same!
            .unwrap()
            .bytes()
            .await;
        // println!("{:?}", response);
        let mut f = std::fs::File::create(path_to_ncbi_dataset)?;
        f.write_all(&response?)?;
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
fn accession_to_lineage(
    df: &DataFrame,
    accession_to_lineage_df: &DataFrame,
    file_path_clade_to_lineages: &String,
) -> Result<DataFrame> {
    dbg!(&df.shape());
    dbg!(&accession_to_lineage_df.shape());

    //read lineage to clade resource
    let lin_to_clade = CsvReader::from_path(file_path_clade_to_lineages)?
        .with_delimiter(b'\t')
        .has_header(true)
        .finish()?;
    dbg!(&lin_to_clade.shape());

    //clone the input df
    let mut renamed_df = df.clone();

    //prepare iterables of columns beforehand
    let mut iter_pangolin = accession_to_lineage_df["Virus Pangolin Classification"]
        .utf8()
        .unwrap()
        .into_iter();
    let mut iter_accession = accession_to_lineage_df["Accession"]
        .utf8()
        .unwrap()
        .into_iter();
    let mut lineages_in_resource = lin_to_clade["lineage"].utf8().unwrap();
    let mut clades_in_resource = lin_to_clade["clade"].utf8().unwrap();

    //rename accessions to clades, using pango lineages.
    iter_pangolin
        .into_iter()
        .zip(iter_accession)
        .for_each(|(png, acc)| {
            let mut rename_str = format!("no clade");
            lineages_in_resource
                .into_iter()
                .zip(clades_in_resource.into_iter())
                .for_each(|(lng, cld)| {
                    if png == lng {
                        rename_str = format!("{}", cld.unwrap());
                    }
                });
            renamed_df.rename(acc.unwrap(), &rename_str);
        });
    dbg!(&renamed_df.shape());
    Ok(renamed_df)
}

fn download_accession(accession: String, path: &PathBuf) -> Result<ExitStatus> {
    let download_genomes_status = {
        Command::new("datasets")
            .arg("download")
            .arg("virus")
            .arg("genome")
            .arg("accession")
            .arg(&accession)
            .arg("--filename")
            .arg(&path)
            .status()
            .expect("failed to execute the download process")
    };

    Ok(download_genomes_status)
}
