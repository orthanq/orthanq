use anyhow::Result;
use derive_builder::Builder;
use std::env::temp_dir;
use std::path::PathBuf;
use std::process::Command;
use tempfile::NamedTempFile;

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
        let index = {
            Command::new("bwa")
                .arg("index")
                .arg("-p")
                .arg(&temp_index.path())
                .arg("-a")
                .arg("bwtsw") //-a bwtsw' does not work for short genomes, lineage quantification?
                .arg(self.genome.clone())
                .status()
                .expect("failed to execute indexing process")
        };
        println!("The index was created successfully: {}", index);

        //perform the alignment for paired end reads
        let mut temp_aligned = NamedTempFile::new()?;

        //find sample name of one of the fastq files from the read pair
        let sample_name = &self.reads[0].file_stem().unwrap();

        //insert read_group info from the sample names
        let read_group = format!("@RG\tID:{:?}\tSM:{:?}", sample_name, sample_name);

        //align reads to the bwa index
        let align = {
            Command::new("bwa")
                .arg("mem")
                .arg("-t")
                .arg("20")
                .arg("-R")
                .arg(&read_group)
                .arg(&temp_index.path())
                .arg(&self.reads[0])
                .arg(&self.reads[1])
                .arg(">")
                .arg(temp_aligned.path())
                .status()
                .expect("failed to execute the alignment process")
        };
        println!("The alignment was done successfully: {}", align);

        //sort the aligned reads by coordinate
        let temp_aligned_sorted = NamedTempFile::new()?;
        let sort = {
            Command::new("samtools")
                .arg("sort")
                .arg(temp_aligned.path())
                .arg("-o")
                .arg(temp_aligned_sorted.path())
                .status()
                .expect("failed to execute the sorting process")
        };
        println!("The sorting was done successfully: {}", align);

        //close the file handle of the named temporary files
        temp_index.close()?;
        temp_aligned.close()?;
        temp_aligned_sorted.close()?;

        Ok(())
    }
}
