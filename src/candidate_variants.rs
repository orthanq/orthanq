use anyhow::Result;
use bio_types::genome::AbstractInterval;
use derive_builder::Builder;
use ndarray::Array2;
use polars::df;
use polars::frame::DataFrame;
use polars::prelude::CsvWriter;
use polars::prelude::NamedFrom;
use polars::prelude::SerWriter;
use polars::series::Series;
use rust_htslib::bcf::header::Header;
use rust_htslib::bcf::{Format, Writer};
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::record::Cigar, bam::Read, faidx};
use std::collections::BTreeMap;
use std::convert::TryInto;
use std::fs;
use std::fs::File;
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
        //self.alignment();
        //TODO: generation of candidate variants from the alignment.

        //1) first read the hla alleles for the locus and the reference genome.
        //let mut sam = bam::Reader::from_path(&"alignment_sorted.sam").unwrap();
        let mut sam = bam::Reader::from_path(&"first10.bam").unwrap();
        let reference_genome = faidx::Reader::from_path(&self.genome).unwrap();

        //2) loop over the alignment file to record the snv and small indels.
        let mut candidate_variants: BTreeMap<
            (String, usize, String, String, &str, &str, &str, &str),
            Vec<usize>,
        > = BTreeMap::new();
        let mut seq_names: Vec<String> = Vec::new();
        let mut locations: Vec<(usize, String, usize, usize)> = Vec::new();
        let mut j: usize = 0; //the iterator that stores the index of name of the allele
        for seq in sam.records() {
            let seq = seq?;
            if !seq.is_secondary() {
                seq_names.push(String::from_utf8(seq.qname().to_vec()).unwrap());
                locations.push((
                    j,
                    seq.contig().to_string(),
                    (seq.reference_start() + 1).try_into().unwrap(),
                    (seq.reference_end() + 1).try_into().unwrap(),
                ));
                //store haplotype locations on the aligned genome
                let mut rcount: i64 = 0; //count for reference position
                let mut scount: i64 = 0; //count for mapped sequence position
                let mut query_alignment_start = 0;
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
                                let spos = scount + query_alignment_start + i;

                                let chrom = seq.contig().to_string();
                                let pos = rpos;
                                let ref_base = reference_genome
                                    .fetch_seq_string(&seq.contig().to_string(), rpos - 1, rpos - 1)
                                    .unwrap();
                                let alt_base = (seq.seq()[spos as usize] as char).to_string();
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
                            let spos = scount + query_alignment_start - 1; //the last matching or mismatch position

                            let chrom = seq.contig().to_string();
                            let pos = rpos;
                            let ref_base = reference_genome
                                .fetch_seq_string(&seq.contig().to_string(), rpos - 1, rpos - 1)
                                .unwrap();
                            let alt_sequence = Vec::from_iter(spos..spos + num + 1)
                                .iter()
                                .map(|pos| seq.seq()[*pos as usize] as char)
                                .collect::<String>();
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
                            let spos = scount + query_alignment_start - 1; //the last matching or mismatch position

                            let chrom = seq.contig().to_string();
                            let pos = rpos;
                            let ref_sequence = reference_genome
                                .fetch_seq_string(
                                    &seq.contig().to_string(),
                                    rpos - 1,
                                    rpos + (num as usize) - 1,
                                )
                                .unwrap();
                            let alt_base = (seq.seq()[spos as usize] as char).to_string();
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
                            query_alignment_start += num * 2;
                            scount -= num;
                        }
                        _ => (),
                    }
                }
                j += 1;
            }
        }

        //construct the first array having the same number of rows and columns as candidate variants map and locate the genotypes for each haplotype.
        let mut genotypes_array = Array2::<i64>::zeros((candidate_variants.len(), j));
        candidate_variants
            .iter()
            .enumerate()
            .for_each(|(i, (_, haplotype_indices))| {
                haplotype_indices
                    .iter()
                    .for_each(|haplotype| genotypes_array[[i, *haplotype]] = 1)
            });
        dbg!(&genotypes_array.shape());
        //construct the second array and locate the loci information of variants for each haplotype
        let mut loci_array = Array2::<i64>::zeros((candidate_variants.len(), j));
        candidate_variants.iter().enumerate().for_each(
            |(i, ((chrom_candidates, pos, _, _, _, _, _, _), _))| {
                locations
                    .iter()
                    .for_each(|(haplotype_index, chrom_locations, start, end)| {
                        if chrom_candidates == chrom_locations && (pos >= start && pos <= end) {
                            loci_array[[i, *haplotype_index]] = 1;
                        }
                    });
            },
        );

        //convert genotype and loci arrays to dataframe.
        //first initialize the dataframes with index columns which contain variant information.
        let mut genotype_df: DataFrame = df!("Index" => candidate_variants
        .iter()
        .map(|((chrom, pos, ref_base, alt_base, n1, n2, n3, gt_c), _)| {
            format!(
                "{}, {}, {}, {}, {}, {}, {}, {}",
                chrom, pos, ref_base, alt_base, n1, n2, n3, gt_c
            )
        })
        .collect::<Series>())?;
        let mut loci_df = genotype_df.clone();
        //second, fill the genotypes dataframe.
        for (index, haplotype_name) in seq_names.iter().enumerate() {
            let query = genotype_df.select_series(&[haplotype_name.as_str()]);
            //the following has to be done in case the haplotype name already exists in the dataframe because
            //of the supplementary alignment records and it needs to be checked to prevent overwriting of the columns.
            //if the haplotype already exists, then the new variants coming from the supplementary alignments
            //are added to the haplotype column to keep all under the same haplotype.
            match query {
                Ok(answer) => genotype_df.with_column(
                    &Series::new(
                        &haplotype_name.as_str(),
                        genotypes_array.column(index).to_vec(),
                    ) + answer.iter().nth(0).unwrap(),
                )?,
                Err(_) => genotype_df.with_column(Series::new(
                    &haplotype_name.as_str(),
                    genotypes_array.column(index).to_vec(),
                ))?,
            };
            //the above operation has to be done for the loci dataframe as well.
            let query = loci_df.select_series(&[haplotype_name.as_str()]);
            match query {
                Ok(answer) => loci_df.with_column(
                    &Series::new(&haplotype_name.as_str(), loci_array.column(index).to_vec())
                        + answer.iter().nth(0).unwrap(),
                )?,
                Err(_) => loci_df.with_column(Series::new(
                    &haplotype_name.as_str(),
                    loci_array.column(index).to_vec(),
                ))?,
            };
        }
        //todo: concatenate the genotype and loci dfs, sort, arrange columns, and write to vcf.

        // let mut output_file: File = File::create("out.csv").unwrap();
        // CsvWriter::new(&mut output_file)
        //     .has_header(false)
        //     .finish(&mut dnm)
        //     .unwrap();

        // // Create minimal VCF header with a single sample
        // let mut header = Header::new();
        // header.push_sample("sample".as_bytes());

        // // Write uncompressed VCF to stdout with above header and get an empty record
        // let mut header = Header::new();
        // let header_contig_line = r#"##contig=<ID=1>"#;
        // header.push_record(header_contig_line.as_bytes());
        // let header_contig_line = r#"##contig=<ID=6>"#;
        // header.push_record(header_contig_line.as_bytes());
        // let header_contig_line = r#"##contig=<ID=8>"#;
        // header.push_record(header_contig_line.as_bytes());
        // let header_contig_line = r#"##contig=<ID=9>"#;
        // header.push_record(header_contig_line.as_bytes());
        // let header_contig_line = r#"##contig=<ID=11>"#;
        // header.push_record(header_contig_line.as_bytes());
        // let header_contig_line = r#"##contig=<ID=16>"#;
        // header.push_record(header_contig_line.as_bytes());
        // let header_contig_line = r#"##contig=<ID=X>"#;
        // header.push_record(header_contig_line.as_bytes());

        // let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Variant is present in the haplotype (1) or not (0).">"#;
        // header.push_record(header_gt_line.as_bytes());
        // genotypes_df.get_column_names().iter().for_each(|sample_name| header.push_sample(sample_name.as_bytes()));
        // let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();

        // //let mut record = vcf.empty_record();
        // for row in genotype_df.rows() {
        //     let mut record = vcf.empty_record();
        //     // Set chrom and pos to 1 and 7, respectively - note the 0-based positions
        //     let rid = vcf.header().name2rid(row).unwrap();
        //     record.set_rid(Some(rid));
        //     record.set_pos(6);

        //     // Set record genotype to 0|1 - note first allele is always unphased
        //     let alleles = &[GenotypeAllele::Unphased(0), GenotypeAllele::Phased(1)];
        //     record.push_genotypes(alleles).unwrap();

        //     // Write record
        //     vcf.write(&record).unwrap()
        // }
        // Writer::write(&mut vcf, &record);
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
