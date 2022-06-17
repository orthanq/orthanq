use anyhow::Result;
use bio_types::genome::AbstractInterval;
use bio_types::sequence::SequenceRead;
use derive_builder::Builder;
use rust_htslib::{
    bam,
    bam::ext::BamRecordExtensions,
    bam::record::{Cigar, CigarString},
    bam::Read,
    faidx,
};
use std::collections::HashMap;
use std::convert::TryInto;
use std::fs;
use std::iter::FromIterator;
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
        let mut sam = bam::Reader::from_path(&"alignment_sorted.sam").unwrap();
        let mut reference_genome = faidx::Reader::from_path(&self.genome).unwrap();

        //2) loop over the alignment file to record the snv and small indels.
        let mut candidate_variants: HashMap<
            (String, usize, String, String, &str, &str, &str, &str),
            Vec<i64>,
        > = HashMap::new();
        let mut seq_names: Vec<String> = Vec::new();
        let mut locations: Vec<(i64, String, i64, i64)> = Vec::new();
        let mut j: i64 = 0; //the iterator that stores the index of name of the allele

        for seq in sam.records() {
            let seq = seq?;
            if !seq.is_secondary() {
                seq_names.push(String::from_utf8(seq.qname().to_vec()).unwrap());
                locations.push((
                    j,
                    seq.contig().to_string(),
                    seq.reference_start() + 1,
                    seq.reference_end() + 1,
                ));
                //store haplotype locations on the aligned genome
                let mut rcount: i64 = 0; //count for reference position
                let mut scount: i64 = 0; //count for mapped sequence position

                for cigar_view in &seq.cigar() {
                    //store the detected variants for each record.
                    match cigar_view {
                        Cigar::Equal(num) => {
                            let num = i64::from(*num);
                            rcount += num;
                            scount += num;
                        }
                        Cigar::Diff(num) => {
                            let num = i64::from(*num);
                            for i in 0..num {
                                let rpos = (rcount + seq.reference_start() + 1 + i) as usize;
                                let spos = scount + i;

                                let chrom = seq.contig().to_string();
                                let pos = rpos;
                                let ref_base = reference_genome
                                    .fetch_seq_string(&seq.contig().to_string(), rpos - 1, rpos)
                                    .unwrap();
                                dbg!(&ref_base);
                                dbg!(&spos);
                                let alt_base = (seq.seq()[spos as usize] as char).to_string();
                                dbg!(&alt_base);
                                candidate_variants
                                    .entry((
                                        chrom.clone(),
                                        pos,
                                        ref_base.clone(),
                                        alt_base.clone(),
                                        ".",
                                        ".",
                                        ".",
                                        "GT:C",
                                    ))
                                    .or_insert(vec![]);
                                let mut haplotypes = candidate_variants
                                    .get(&(
                                        chrom.clone(),
                                        pos,
                                        ref_base.clone(),
                                        alt_base.clone(),
                                        ".",
                                        ".",
                                        ".",
                                        "GT:C",
                                    ))
                                    .unwrap()
                                    .clone();
                                haplotypes.push(j); //appending j informs the dict about the allele names eventually
                                candidate_variants.insert(
                                    (chrom, pos, ref_base, alt_base, ".", ".", ".", "GT:C"),
                                    haplotypes,
                                );
                            }
                            rcount += num; //to add mismatch length to the count
                            scount += num;
                        }
                        Cigar::Ins(num) => {
                            let num = i64::from(*num);

                            let rpos = (rcount + seq.reference_start()) as usize; //the last matching or mismatch position
                            let spos = scount - 1; //the last matching or mismatch position

                            let chrom = seq.contig().to_string();
                            let pos = rpos;
                            let ref_base = reference_genome
                                .fetch_seq_string(&seq.contig().to_string(), rpos - 1, rpos)
                                .unwrap();
                            dbg!(&ref_base);
                            let alt_sequence = Vec::from_iter(spos..spos + num + 1)
                                .iter()
                                .map(|pos| seq.seq()[*pos as usize] as char)
                                .collect::<String>();
                            dbg!(&alt_sequence);
                            candidate_variants
                                .entry((
                                    chrom.clone(),
                                    pos,
                                    ref_base.clone(),
                                    alt_sequence.clone(),
                                    ".",
                                    ".",
                                    ".",
                                    "GT:C",
                                ))
                                .or_insert(vec![]);
                            let mut haplotypes = candidate_variants
                                .get(&(
                                    chrom.clone(),
                                    pos,
                                    ref_base.clone(),
                                    alt_sequence.clone(),
                                    ".",
                                    ".",
                                    ".",
                                    "GT:C",
                                ))
                                .unwrap()
                                .clone();
                            haplotypes.push(j); //appending j informs the dict about the allele names eventually
                            candidate_variants.insert(
                                (chrom, pos, ref_base, alt_sequence, ".", ".", ".", "GT:C"),
                                haplotypes,
                            );
                            rcount; // no change
                            scount += num;
                        }
                        Cigar::Del(num) => {
                            let num = i64::from(*num);

                            let rpos = (rcount + seq.reference_start()) as usize; //the last matching or mismatch position
                            let spos = scount - 1; //the last matching or mismatch position

                            let chrom = seq.contig().to_string();
                            let pos = rpos;
                            dbg!(&pos);
                            let ref_sequence = reference_genome
                                .fetch_seq_string(
                                    &seq.contig().to_string(),
                                    rpos - 1,
                                    rpos + (num as usize),
                                )
                                .unwrap();
                            dbg!(&ref_sequence);
                            let alt_base = (seq.seq()[spos as usize] as char).to_string();
                            dbg!(&alt_base);
                            candidate_variants
                                .entry((
                                    chrom.clone(),
                                    pos,
                                    ref_sequence.clone(),
                                    alt_base.clone(),
                                    ".",
                                    ".",
                                    ".",
                                    "GT:C",
                                ))
                                .or_insert(vec![]);
                            let mut haplotypes = candidate_variants
                                .get(&(
                                    chrom.clone(),
                                    pos,
                                    ref_sequence.clone(),
                                    alt_base.clone(),
                                    ".",
                                    ".",
                                    ".",
                                    "GT:C",
                                ))
                                .unwrap()
                                .clone();
                            haplotypes.push(j); //appending j informs the dict about the allele names eventually
                            candidate_variants.insert(
                                (chrom, pos, ref_sequence, alt_base, ".", ".", ".", "GT:C"),
                                haplotypes,
                            );

                            rcount += num;
                            scount;
                        }
                        Cigar::SoftClip(num) => {
                            let num = i64::from(*num);
                            scount -= num;
                        }
                        _ => (),
                    }
                    j += 1;
                }
            }
        }
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
