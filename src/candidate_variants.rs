use anyhow::Result;
use bio_types::genome::AbstractInterval;
use derive_builder::Builder;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::Read, faidx};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

#[allow(dead_code)]
#[derive(Builder, Clone)]
pub struct Caller {
    alleles: PathBuf,
    genome: PathBuf,
}
impl Caller {
    pub fn call(&self) -> Result<()> {
        self.alignment();
        //TODO: generation of candidate variants from the alignment.

        //1) first read the hla alleles for the locus and the reference genome.
        let mut sam = bam::Reader::from_path(&"alignment_sorted.sam").unwrap(); //separate alignments for each locus.
        let mut reference_genome = faidx::Reader::from_path(&"DQA1_filtered.fasta"); //normally hs_genome.fasta

        //2) loop over the alignment file to record the snv and small indels.
        let mut candidate_variants: HashMap<(), Vec<usize>> = HashMap::new();
        let mut seq_names: Vec<String> = Vec::new();
        let mut locations: Vec<(&i64, String, i64, i64)> = Vec::new();
        let mut j: i64 = 0; //the iterator that stores the index of name of the allele

        for seq in sam.records() {
            let seq = seq?;
            if !seq.is_secondary() {
                seq_names.push(String::from_utf8(seq.qname().to_vec()).unwrap());
                locations.push((
                    &j,
                    seq.contig().to_string(),
                    seq.reference_start() + 1,
                    seq.reference_end() + 1,
                )); //store haplotype locations on the aligned genome
                let rcount: i64 = 0; //count for reference position
                let scount: i64 = 0; //count for mapped sequence position
                for cigar_view in seq.cigar() {
                    //store the detected variants for each record.
                    todo!();
                }
            }
        }
        // dbg!(&seq_names);
        // dbg!(&locations);
        Ok(())
    }
    #[allow(dead_code)]
    fn alignment(&self) -> Result<()> {
        let genome_name = format!(
            "{}{}",
            self.genome.file_name().unwrap().to_str().unwrap(),
            ".mmi"
        );
        let _alleles_name = format!("{}", self.alleles.file_name().unwrap().to_str().unwrap());
        let index = {
            Command::new("minimap2")
                .arg("-d")
                .arg(&genome_name)
                .arg(self.genome.clone())
                .status()
                .expect("failed to execute indexing process")
        };
        println!("indexing process finished with: {}", index);

        let align = {
            Command::new("minimap2")
                .args(["-a", "-t", "36"])
                .arg(&genome_name)
                .arg(self.alleles.clone())
                .output()
                .expect("failed to execute alignment process")
        };
        let stdout = String::from_utf8(align.stdout).unwrap();
        println!("alignment process finished!");
        fs::write("alignment.sam", stdout).expect("Unable to write file");

        //sort and convert the sam to bam
        let sort = {
            Command::new("samtools")
                .arg("sort")
                .arg("alignment.sam")
                .output()
                .expect("failed to execute alignment process")
        };
        let stdout = sort.stdout;
        println!("sorting process finished!");
        fs::write("alignment_sorted.sam", stdout).expect("Unable to write file");
        Ok(())
    }
}
