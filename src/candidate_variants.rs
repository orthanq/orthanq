use anyhow::Result;
use bio::io::fasta::Reader;
use bio_types::genome::AbstractInterval;
use csv::Reader as CsvReader;
use derive_builder::Builder;
use ndarray::Array2;
use ordered_float::NotNan;
use polars::{
    datatypes::Utf8Chunked,
    df,
    frame::DataFrame,
    prelude::{CsvWriter, IntoSeries, NamedFrom, SerWriter},
    series::Series,
};
use quick_xml::events::Event;
use quick_xml::reader::Reader as xml_reader;
use rust_htslib::bcf::{header::Header, record::GenotypeAllele, Format, Writer};
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::record::Cigar, bam::Read, faidx};
use serde::Deserialize;
use std::collections::{BTreeMap, HashMap};
use std::convert::TryInto;
use std::error::Error;
use std::io;
use std::iter::FromIterator;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;
use std::{fs, fs::File};

#[allow(dead_code)]
#[derive(Builder, Clone)]
pub struct Caller {
    alleles: PathBuf,
    genome: PathBuf,
    xml: PathBuf,
    allele_freq: PathBuf,
    wes: bool,
    wgs: bool,
    output: Option<PathBuf>,
}
impl Caller {
    pub fn call(&self) -> Result<()> {
        //prepare the map to look up which alleles are confirmed and unconfirmed and g codes available
        //IMGT/HLA version for xml file is 3.35, set this to the same version in the evaluation workflow
        //and only include alleles that have AF >0.05
        let af_path = "allele_frequencies.csv".to_string();
        let (confirmed_alleles, unconfirmed_alleles) =
            confirmed_alleles(&self.xml, &self.allele_freq).unwrap();

        //write loci to separate fasta files (Confirmed and alleles that have g codes available)
        // self.write_to_fasta(&confirmed_alleles)?;

        //align and sort
        self.alignment(); //copied from exclude folder to include folder so i commented this out for a short while.

        //1) first read the hla alleles for the locus and the reference genome.
        let mut sam = bam::Reader::from_path(&"alignment_sorted.sam").unwrap();
        let reference_genome = faidx::Reader::from_path(&self.genome).unwrap();

        //2) loop over the alignment file to record the snv and small indels.
        let mut candidate_variants: BTreeMap<(String, usize, String, String), Vec<usize>> =
            BTreeMap::new();
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
                                    .entry((chrom.clone(), pos, ref_base.clone(), alt_base.clone()))
                                    .or_insert(vec![]);
                                let mut haplotypes = candidate_variants
                                    .get(&(chrom.clone(), pos, ref_base.clone(), alt_base.clone()))
                                    .unwrap()
                                    .clone();
                                haplotypes.push(j); //appending j informs the dict about the allele names eventually
                                candidate_variants
                                    .insert((chrom, pos, ref_base, alt_base), haplotypes);
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
                                .entry((chrom.clone(), pos, ref_base.clone(), alt_sequence.clone()))
                                .or_insert(vec![]);
                            let mut haplotypes = candidate_variants
                                .get(&(chrom.clone(), pos, ref_base.clone(), alt_sequence.clone()))
                                .unwrap()
                                .clone();
                            haplotypes.push(j); //appending j informs the dict about the allele names eventually
                            candidate_variants
                                .insert((chrom, pos, ref_base, alt_sequence), haplotypes);
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
                                .entry((chrom.clone(), pos, ref_sequence.clone(), alt_base.clone()))
                                .or_insert(vec![]);
                            let mut haplotypes = candidate_variants
                                .get(&(chrom.clone(), pos, ref_sequence.clone(), alt_base.clone()))
                                .unwrap()
                                .clone();
                            haplotypes.push(j); //appending j informs the dict about the allele names eventually
                            candidate_variants
                                .insert((chrom, pos, ref_sequence, alt_base), haplotypes);

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

        //construct the second array and locate the loci information of variants for each haplotype
        let mut loci_array = Array2::<i64>::zeros((candidate_variants.len(), j));
        candidate_variants.iter().enumerate().for_each(
            |(i, ((chrom_candidates, pos, _, _), _))| {
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
        .map(|((chrom, pos, ref_base, alt_base), _)| {
            format!(
                "{},{},{},{}",
                chrom, pos, ref_base, alt_base
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
        // there is no need to sort the map as we used BTreeMap to store the variants.

        //insert an ID column as many as the number of rows, right after the Index column
        genotype_df.insert_at_idx(
            1,
            Series::new("ID", Vec::from_iter(0..genotype_df.shape().0 as i64)),
        )?;
        loci_df.insert_at_idx(
            1,
            Series::new("ID", Vec::from_iter(0..loci_df.shape().0 as i64)),
        )?;

        //Unconfirmed alleles are removed from both dataframes
        unconfirmed_alleles.iter().for_each(|unconf_allele| {
            genotype_df.drop_in_place(unconf_allele);
            loci_df.drop_in_place(unconf_allele);
        });

        //Split up the variants per locus and depending on --wes samples, restrict the columns to protein level.
        //TODO: enable writing for wgs samples.
        if self.wes {
            let genotype_df = self.split_haplotypes(&mut genotype_df)?;
            let loci_df = self.split_haplotypes(&mut loci_df)?;
            //write locus-wise vcf files.
            self.write_loci_to_vcf(&genotype_df, &loci_df)?;
        }
        Ok(())
    }

    fn split_haplotypes(&self, variant_table: &mut DataFrame) -> Result<DataFrame> {
        let hla_alleles_rdr = bio::io::fasta::Reader::from_file(&self.alleles)?;
        let mut allele_digit_table = HashMap::new();
        for record in hla_alleles_rdr.records() {
            let record = record.unwrap();
            let id = record.id();
            let desc = record.desc().unwrap();

            let six_digit = desc.split_whitespace().collect::<Vec<&str>>()[0];
            let splitted = six_digit.split(":").collect::<Vec<&str>>();
            if splitted.len() < 2 {
                //some alleles e.g. MICA may not have the full nomenclature, i.e. 6 digits
                allele_digit_table.insert(id.to_string(), splitted[0].to_string());
            } else if splitted.len() < 3 {
                //some alleles e.g. MICA may not have the full nomenclature, i.e. 6 digits
                allele_digit_table
                    .insert(id.to_string(), format!("{}:{}", splitted[0], splitted[1]));
            } else {
                allele_digit_table.insert(
                    id.to_string(),
                    format!("{}:{}:{}", splitted[0], splitted[1], splitted[2]),
                );
                //first two
            }
            // else {
            //     allele_digit_table.insert(
            //         id.to_string(),
            //         format!(
            //             "{}:{}:{}:{}",
            //             splitted[0], splitted[1], splitted[2], splitted[3]
            //         ),
            //     );
            // }
        }

        let mut new_df = DataFrame::new(vec![
            variant_table["Index"].clone(),
            variant_table["ID"].clone(),
        ])?;
        for (column_index, column_name) in
            variant_table.get_column_names().iter().skip(2).enumerate()
        {
            let protein_level = &allele_digit_table[&column_name.to_string()];
            if new_df.get_column_names().contains(&protein_level.as_str()) {
                let existing_column = &new_df[protein_level.as_str()];
                let updated_column = existing_column + &variant_table[*column_name];
                new_df.with_column(updated_column)?;
            } else {
                new_df.replace_or_add(&protein_level, variant_table[*column_name].clone())?;
            }
        }
        let names: Vec<String> = new_df.get_column_names_owned().iter().cloned().collect();
        for (column_index, column_name) in names.iter().enumerate().skip(2) {
            let corrected_column = new_df[column_index]
                .i64()
                .unwrap()
                .into_iter()
                .map(|num| if num > Some(1) { Some(1) } else { num })
                .collect::<Series>();
            new_df.replace_or_add(column_name, corrected_column)?;
        }
        Ok(new_df)
    }

    #[allow(dead_code)]
    fn alignment(&self) -> Result<()> {
        let genome_name = format!(
            "{}{}",
            self.genome.file_name().unwrap().to_str().unwrap(),
            ".mmi"
        );
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
                .args(["-a", "-t", "36", "--eqx", "--MD"])
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

    fn write_loci_to_vcf(&self, variant_table: &DataFrame, loci_table: &DataFrame) -> Result<()> {
        let names = variant_table
            .get_column_names()
            .iter()
            .cloned()
            .collect::<Vec<&str>>();

        // for locus in vec![
        //     "A", "DPA1", "DRB4", "V", "B", "DPB1", "DRB5", "W", "C", "DQA1", "E", "DQA2", "F", "S",
        //     "DMA", "DQB1", "G", "TAP1", "DMB", "DRA", "HFE", "TAP2", "DOA", "DRB1", "T", "DOB",
        //     "DRB3", "MICA", "U", //all the non-pseudogenes
        // ]
        for locus in vec!["A", "B", "C", "DQA1", "DQB1"] {
            let mut locus_columns = vec!["Index".to_string(), "ID".to_string()];
            for column_name in names.iter().skip(2) {
                let splitted = column_name.split("*").collect::<Vec<&str>>();
                if splitted[0] == locus {
                    // just as an example locus for the start.
                    locus_columns.push(column_name.to_string());
                }
            }

            let variant_table = variant_table.select(&locus_columns)?;
            let loci_table = loci_table.select(&locus_columns)?;

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
            //fs::create_dir_all(self.output.as_ref().unwrap())?;
            let mut vcf = Writer::from_path(
                format!(
                    "{}.vcf",
                    self.output.as_ref().unwrap().join(locus).display()
                ),
                &header,
                true,
                Format::Vcf,
            )
            .unwrap();

            let id_iter = variant_table["ID"].i64().unwrap().into_iter();
            for row_index in 0..variant_table.height() {
                let mut record = vcf.empty_record();
                let mut variant_iter = variant_table["Index"].utf8().unwrap().into_iter();
                let mut id_iter = variant_table["ID"].i64().unwrap().into_iter();
                let splitted = variant_iter
                    .nth(row_index)
                    .unwrap()
                    .unwrap()
                    .split(",")
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
                        .i64()
                        .unwrap()
                        .into_iter()
                        .nth(row_index)
                        .unwrap()
                        .unwrap();
                    let gt = GenotypeAllele::Unphased(gt.try_into().unwrap());
                    all_gt.push(gt);
                }
                record.push_genotypes(&all_gt).unwrap();

                //push loci
                let mut all_c = Vec::new();
                for column_index in 2..loci_table.width() {
                    //it doesnt have the ID column so that it starts from 1
                    let c = loci_table[column_index]
                        .i64()
                        .unwrap()
                        .into_iter()
                        .nth(row_index)
                        .unwrap()
                        .unwrap();
                    all_c.push(c as i32);
                }
                record.push_format_integer(b"C", &all_c)?;
                vcf.write(&record).unwrap();
            }
        }
        Ok(())
    }
    //Write_to_fasta() function collects confirmed allele list from confirmed_alleles() function.
    //Then it generates Fasta files for those alleles for each loci (necessary for quantifications e.g. kallisto, salmon)
    fn write_to_fasta(&self, confirmed_alleles: &Vec<String>) -> Result<()> {
        for locus in vec!["A", "B", "C", "DQA1", "DQB1"] {
            fs::create_dir_all(self.output.as_ref().unwrap())?;
            let file = fs::File::create(format!(
                "{}.fasta",
                self.output.as_ref().unwrap().join(locus).display()
            ))
            .unwrap();
            let handle = io::BufWriter::new(file);
            let mut writer = bio::io::fasta::Writer::new(handle);
            for record in bio::io::fasta::Reader::from_file(&self.alleles)?.records() {
                let record = record.unwrap();
                let id = record.id();
                let description = record.desc().unwrap();
                let splitted = description.split("*").collect::<Vec<&str>>();
                if splitted[0] == locus {
                    if confirmed_alleles.contains(&id.to_string()) {
                        let new_record = bio::io::fasta::Record::with_attrs(
                            description.split_whitespace().collect::<Vec<&str>>()[0].clone(),
                            Some(id.clone()),
                            record.seq(),
                        );
                        writer
                            .write_record(&new_record)
                            .ok()
                            .expect("Error writing record.");
                    }
                }
            }
        }
        Ok(())
    }
}

#[derive(Debug, Deserialize)]
struct Record {
    var: String,
    population: String,
    frequency: NotNan<f64>,
}

//Confirmed_alleles() function finds HLA alleles that are "Confirmed" and "Unconfirmed".
//In the end the unconfirmed vector contains "Unconfirmed" alleles

//update: remove alleles that are not AF <0.05 in at least one population.

fn confirmed_alleles(xml_path: &PathBuf, af_path: &PathBuf) -> Result<(Vec<String>, Vec<String>)> {
    let mut reader = xml_reader::from_file(&xml_path)?;
    reader.trim_text(true);
    let mut buf = Vec::new();
    let mut alleles: Vec<String> = Vec::new();
    let mut allele_names: Vec<String> = Vec::new();
    let mut confirmed: Vec<String> = Vec::new();
    let mut hla_g_groups: HashMap<i32, String> = HashMap::new(); //some hla alleles dont have g groups information in the xml file.
    let mut alleles_indices: Vec<i32> = Vec::new();
    let mut groups_indices: Vec<i32> = Vec::new();
    let mut counter = 0;
    loop {
        match reader.read_event_into(&mut buf) {
            Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
            Ok(Event::Eof) => break,
            Ok(Event::Start(e)) => match e.name().as_ref() {
                b"allele" => {
                    alleles.push(
                        e.attributes()
                            .map(|a| String::from_utf8(a.unwrap().value.to_vec()))
                            .collect::<Vec<_>>()[0]
                            .as_ref()
                            .unwrap()
                            .to_string(),
                    ); //index 0 holds the allele id
                    allele_names.push(
                        e.attributes()
                            .map(|a| String::from_utf8(a.unwrap().value.to_vec()))
                            .collect::<Vec<_>>()[1]
                            .as_ref()
                            .unwrap()
                            .split("-") //"HLA-" don't take the HLA prefix
                            .collect::<Vec<&str>>()[1]
                            .to_string(),
                    ); //index 1 holds the allele name
                    alleles_indices.push(counter.clone());
                    counter += 1;
                }
                _ => (),
            },
            Ok(Event::Empty(e)) => match e.name().as_ref() {
                b"releaseversions" => confirmed.push(
                    e.attributes()
                        .map(|a| String::from_utf8(a.unwrap().value.to_vec()))
                        .collect::<Vec<_>>()[4]
                        .as_ref()
                        .unwrap()
                        .to_string(), //index 4 holds the Confirmed info
                ),
                b"hla_g_group" => {
                    groups_indices.push(counter.clone());
                    hla_g_groups.insert(
                        counter.clone(),
                        e.attributes()
                            .map(|a| String::from_utf8(a.unwrap().value.to_vec()))
                            .collect::<Vec<_>>()[0]
                            .as_ref()
                            .unwrap()
                            .to_string(), //index 0 holds the status info
                    );
                }
                _ => (),
            },
            _ => (),
        }
        // if we don't keep a borrow elsewhere, we can clear the buffer to keep memory usage low
        buf.clear();
    }
    assert_eq!(alleles.len(), confirmed.len());
    //create a map for allele ids and allele names
    let mut unconfirmed_alleles = alleles
        .iter()
        .zip(allele_names.iter())
        .zip(confirmed.iter())
        .filter(|((id, name), c)| c == &"Unconfirmed")
        .map(|((id, name), c)| (format!("HLA:{}", id.clone()), name.clone()))
        .collect::<HashMap<String, String>>();
    //write confirmed alleles to file
    let mut confirmed_alleles = alleles
        .iter()
        .zip(allele_names.iter())
        .filter(|(id, name)| !unconfirmed_alleles.contains_key(&format!("HLA:{}", id)))
        .map(|(id, name)| (format!("HLA:{}", id.clone()), name.clone()))
        .collect::<HashMap<String, String>>();
    //include only alleles that have >0.05 AF in at least one population
    let mut unconfirmed_alleles = unconfirmed_alleles.keys().cloned().collect::<Vec<String>>();
    dbg!(&unconfirmed_alleles.len());
    dbg!(&confirmed_alleles.len());
    let mut confirmed_alleles_clone = confirmed_alleles.clone();
    let mut allele_freq_rdr = CsvReader::from_path(af_path)?;
    let mut to_be_included = Vec::new();
    dbg!(confirmed_alleles_clone.len());
    confirmed_alleles_clone.retain(|&_, y| {
        y.starts_with("A")
            || y.starts_with("B")
            || y.starts_with("C")
            || y.starts_with("DQA1")
            || y.starts_with("DQB1")
    });
    dbg!(confirmed_alleles_clone.len());

    allele_freq_rdr.deserialize().for_each(|result| {
        let record: Record = result.unwrap();
        confirmed_alleles_clone.iter().for_each(|(id, name)| {
            let splitted = name.split(":").collect::<Vec<&str>>(); //DQB1*05:02:01 is alone not an allele name, but DQB1*05:02:01:01, DQB1*05:02:01:02.. are.
            let first_two = format!("{}:{}", splitted[0], splitted[1]);
            //B*39:06:01 is below 0.05 but 39:06 not and this allele is one of the true genotypes of a sample in the ground truths so we should include following lines. a direct match is not preferred by the authors in the ground truth.
            //they include alleles that do not have a direct name match, rather the ones starting with the first two fields.
            if &record.var == name {
                if record.frequency > NotNan::new(0.05).unwrap() {
                    to_be_included.push(id.clone());
                }
            } else if &record.var == &first_two && record.frequency > NotNan::new(0.05).unwrap() {
                to_be_included.push(id.clone());
            }
        });
    });
    dbg!(&to_be_included.len());
    dbg!(&to_be_included);
    let below_criterium: Vec<String> = confirmed_alleles_clone
        .iter()
        .filter(|(id, _)| !to_be_included.contains(id))
        .map(|(id, _)| id.clone())
        .collect();
    dbg!(&below_criterium.len());
    dbg!(&below_criterium);
    //todo: confirmed_alleles
    let confirmed_alleles = confirmed_alleles.keys().cloned().collect::<Vec<String>>();
    unconfirmed_alleles.extend(below_criterium);
    dbg!(&unconfirmed_alleles.len());
    dbg!(&confirmed_alleles.len());
    Ok((confirmed_alleles, unconfirmed_alleles))
}
