use anyhow::Result;
use derive_builder::Builder;
use std::env::temp_dir;
use std::path::PathBuf;
use std::process::Command;
use tempfile::tempdir;
use tempfile::{NamedTempFile, TempDir};

#[derive(Builder, Clone)]
pub struct Caller {
    genome: PathBuf,
    reads: Vec<PathBuf>,
    output: Option<PathBuf>,
}

impl Caller {
    pub fn call(&self) -> Result<()> {
        //1-) index the linear genome
        //bwa index -p hs_genome -a bwtsw hs_genome.fasta

        //create a temporary file for bwa index and execute bwa index
        let mut temp_index = NamedTempFile::new()?;
        // let index = {
        //     Command::new("bwa")
        //         .arg("index")
        //         .arg("-p")
        //         .arg(&temp_index.path())
        //         .arg("-a")
        //         .arg("bwtsw") //-a bwtsw' does not work for short genomes, lineage quantification?
        //         .arg(self.genome.clone())
        //         .status()
        //         .expect("failed to execute indexing process")
        // };
        // println!("The index was created successfully: {}", index);

        //perform the alignment for paired end reads
        // let mut temp_aligned = NamedTempFile::new()?;

        // Create a directory inside of `std::env::temp_dir()`
        let temp_dir = tempdir()?;

        //find sample name of one of the fastq files from the read pair
        let sample_name = &self.reads[0].file_stem().unwrap().to_str().unwrap();
        //get rid of any underscore (PE reads contain them)
        let sample_name_split = sample_name.split("_").collect::<Vec<&str>>();
        let sample_name = sample_name_split[0];

        //create the output file name in temp directory
        let file_aligned = temp_dir.path().join(format!("{}.bam", sample_name.clone()));
        println!("{}", file_aligned.display());

        //insert read_group info from the sample names
        let read_group = format!(
            "@RG\\tID:{}\\tSM:{}",
            sample_name.clone(),
            sample_name.clone()
        );

        //test the subcommand, for this skip the indexing part
        let index = "/projects/koesterlab/orthanq/orthanq-evaluation/results/bwa-index/hs_genome";

        //align reads to the bwa index
        let align = {
            Command::new("bwa")
                .arg("mem")
                .arg("-t")
                .arg("10")
                .arg("-R")
                .arg(&read_group)
                // .arg(&temp_index.path())
                .arg(index)
                .arg(&self.reads[0])
                .arg(&self.reads[1])
                .arg("-o")
                .arg(&file_aligned)
                // .arg("2>")
                // .arg("log.txt")
                .status()
                .expect("failed to execute the alignment process")
        };
        println!("The alignment was exited with: {}", align);
        println!("{}", file_aligned.display());
        //sort the aligned reads by coordinate

        //create the output file name in temp directory
        let file_aligned_sorted = temp_dir
            .path()
            .join(format!("{}_sorted.bam", sample_name.clone()));

        let sort = {
            Command::new("samtools")
                .arg("sort")
                .arg(file_aligned)
                .arg("-o")
                .arg(&file_aligned_sorted)
                .arg("-@")
                .arg("10")
                .arg("--write-index")
                .status()
                .expect("failed to execute the sorting process")
        };
        println!("The sorting was exited with: {}", align);
        println!("{}", file_aligned_sorted.display());

        //extract reads that map to HLA genes (classical and nonclassical class of genes)

        //create the output file name in temp directory
        let file_extracted = temp_dir
            .path()
            .join(format!("{}_extracted.bam", sample_name.clone()));

        let regions = "resources/regions.bed";

        let extract = {
            Command::new("samtools")
                .arg("view")
                .arg(file_aligned_sorted)
                .arg("-R")
                .arg(regions)
                .arg("--write-index") //??
                .arg("-o")
                .arg(&file_extracted)
                .status()
                .expect("failed to execute the extracting process")
        };
        println!("The extraction was exited with: {}", extract);

        //convert the alignment file to fq

        //create the output file name in temp directory
        let temp_extracted_fq_1 = temp_dir
            .path()
            .join(format!("{}_1.fastq", sample_name.clone()));
        let temp_extracted_fq_2 = temp_dir
            .path()
            .join(format!("{}_2.fastq", sample_name.clone()));

        let bam_to_fq = {
            Command::new("samtools")
                .arg("fastq")
                .arg(file_extracted)
                .arg("-n") //-n for fastq
                .arg("-1")
                .arg(temp_extracted_fq_1)
                .arg("-2")
                .arg(temp_extracted_fq_2)
                .status()
                .expect("failed to execute the extracting process")
        };
        println!("Conversion from BAM to fq was exited with: {}", bam_to_fq);

        //close the file handle of the named temporary files
        temp_index.close()?;
        temp_dir.close()?;

        Ok(())
    }
}
