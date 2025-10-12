use anyhow::Result;

use derive_builder::Builder;

use polars::frame::DataFrame;
use polars::prelude::*;

use crate::candidates::hla::{self};
use crate::candidates::virus::generic::write_to_vcf;

use bio::io::fasta;

use ndarray::Array2;
use petgraph::dot::{Config, Dot};
use petgraph::graph::{Graph, NodeIndex};
use petgraph::prelude::Dfs;

use rust_htslib::faidx;
use seq_io::fasta::Record as OtherRecord;
use serde::Deserialize;

use std::collections::BTreeMap;
use std::collections::HashMap;

use std::fs;
use std::fs::File;
use std::io;

use std::io::Write;
use std::iter::FromIterator;
use std::path::PathBuf;

use std::process::Command;

#[allow(dead_code)]
#[derive(Builder, Clone)]
pub struct Caller {
    output: PathBuf,
    threads: String,
    output_bcf: bool,
}
impl Caller {
    pub fn call(&mut self) -> Result<()> {
        //prepare genome and alleles for the virus (currently, only sars-cov2)

        //first, create the output dir for database setup
        let outdir = &self.output;
        fs::create_dir_all(&outdir)?;

        //download required resources
        let (reference_path, clades_path) = download_resources(&outdir).unwrap();

        //read the clades
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(true)
            .has_headers(true)
            .from_path(clades_path)?;

        //generate mutations per clade and create a graph for clade hierarchy
        let mut clade_mutations: HashMap<String, Vec<String>> = HashMap::new();
        let mut clade_hierarchy: Graph<String, String> = Graph::new();
        for result in rdr.deserialize() {
            let record: Record = result?;
            // println!("{:?}", record);
            let clade = record.clade.clone();
            if record.gene == "nuc" {
                if let Some(alt) = record.alt.clone() {
                    let mutation = format!("{}_{}", record.site, alt);
                    if clade_mutations.contains_key(&clade) {
                        let mut existing_muts = clade_mutations[&clade].clone();
                        existing_muts.push(mutation);
                        clade_mutations.insert(clade.clone(), existing_muts);
                    } else {
                        clade_mutations.insert(clade.clone(), vec![mutation]);
                    }
                }
            } else if record.gene == "clade" {
                //add parent and child nodes if they are not present
                let mut parent_presence = false;
                let mut child_presence = false;

                for idx in 0..clade_hierarchy.node_count() {
                    if clade_hierarchy[NodeIndex::new(idx)] == clade {
                        parent_presence = true
                    }
                    if clade_hierarchy[NodeIndex::new(idx)] == record.site {
                        child_presence = true
                    }
                }
                if parent_presence == false {
                    clade_hierarchy.add_node(clade.clone());
                }
                if child_presence == false {
                    clade_hierarchy.add_node(record.site.clone());
                }

                //add edges as variant sets in direction of parent to child
                let index_parent = clade_hierarchy
                    .node_indices()
                    .find(|i| clade_hierarchy[*i] == clade)
                    .unwrap();
                let index_child = clade_hierarchy
                    .node_indices()
                    .find(|i| clade_hierarchy[*i] == record.site)
                    .unwrap();
                // let mut mutation = format!("{}_{}", record.site, alt);
                clade_hierarchy.add_edge(index_parent, index_child, "".to_string());
            }
        }
        dbg!(&clade_hierarchy);

        //update clade mutations with parent mutations
        let mut clade_mutations_with_parents = clade_mutations.clone();
        clade_mutations.iter().for_each(|(clade, _mutations)| {
            if let Some(index) = &clade_hierarchy
                .node_indices()
                .find(|i| &clade_hierarchy[*i] == clade)
            {
                //iterate over all the nodes that come every entry in clade_mutations
                let mut dfs = Dfs::new(&clade_hierarchy, *index);
                let mut mutations_with_parent = clade_mutations[clade].clone();
                while let Some(visited) = dfs.next(&clade_hierarchy) {
                    let visited = clade_hierarchy[visited].clone();
                    dbg!(&visited);
                    //the first visit is the node itself
                    //insert mutations coming from the parent nodes
                    if visited != *clade && visited != "19A" {
                        //19A is the reference itself - no mutations
                        mutations_with_parent.extend(clade_mutations[&visited].clone());
                    }
                }
                clade_mutations_with_parents.insert(clade.clone(), mutations_with_parent);
            }
        });
        dbg!(&clade_mutations);
        dbg!(&clade_mutations_with_parents);

        //add 19B manually to the graph for visualization purposes (19B -> 19A info not present in the clades.tsv)
        let nineteenb_i = clade_hierarchy.add_node("19B".to_string());
        let nineteena_i = clade_hierarchy
            .node_indices()
            .find(|i| &clade_hierarchy[*i] == "19A")
            .unwrap();
        clade_hierarchy.add_edge(nineteenb_i, nineteena_i, "".to_string());

        //export graph as output
        let graph_dir = outdir.join("sarscov2_lineages.dot");
        let mut f = File::create(graph_dir).unwrap();
        let output = format!(
            "{:?}",
            Dot::with_config(&clade_hierarchy, &[Config::EdgeNoLabel])
        );
        f.write_all(&output.as_bytes())
            .expect("could not write file");

        //create a map and an array as a template of upcoming vcf data structure
        let mut candidate_variants: BTreeMap<(String, usize, String, String), Vec<usize>> =
            BTreeMap::new();
        let reference_genome = faidx::Reader::from_path(&reference_path).unwrap();

        //name of the sarscov2 accession to use as chrom
        let chrom = reference_genome.seq_name(0).unwrap().to_string();

        let mut j: usize = 0; //the iterator that stores the index of clade
                              //some nuc conversions in the clades.tsv has are conversions that support the ref.
                              //examples include: 1) 19A and the variant,	nuc	8782 C, which supports reference (pos 8782 also has C in that position)
                              //2) 22F nuc 15461 A. this could be because 22F is a recombinant of two strains and one of them might have the same mutation 15461

        //remove 22F for now. TODO: update this the following line according to the answers from the maintainers
        clade_mutations_with_parents.remove("22F");

        //loop over the graph to construct candidate variants map
        clade_mutations_with_parents
            .iter()
            .for_each(|(_clade, mutations)| {
                //find record components
                let mut pos: usize = 0;
                let mut alt_base: String = "".to_string();
                mutations.iter().for_each(|m| {
                    dbg!(&m);
                    let splitted = m.split('_').collect::<Vec<&str>>();
                    pos = splitted[0].parse::<usize>().unwrap();
                    alt_base = splitted[1].to_string();

                    let ref_base = reference_genome
                        .fetch_seq_string(&chrom, pos - 1, pos - 1)
                        .unwrap();
                    dbg!(&pos, &ref_base, &alt_base);

                    //update the record in candidate variants if ref base and alt bases are different (update according to the response from Cornelius Roemer)
                    //19A is still being represented because of the count j+=1 in the end and the array was consisting of zeros for all haplotypes.
                    if ref_base != alt_base {
                        candidate_variants
                            .entry((chrom.clone(), pos, ref_base.clone(), alt_base.clone()))
                            .or_insert(vec![]);
                        let mut haplotypes = candidate_variants
                            .get(&(chrom.clone(), pos, ref_base.clone(), alt_base.clone()))
                            .unwrap()
                            .clone();
                        haplotypes.push(j); //appending j informs the dict about the clade names eventually
                        candidate_variants
                            .insert((chrom.clone(), pos, ref_base, alt_base.clone()), haplotypes);
                    }
                });
                j += 1;
            });
        dbg!(&candidate_variants);

        //convert candidate variants to array
        let genotypes_array =
            hla::convert_candidate_variants_to_array(candidate_variants.clone(), j).unwrap();

        //create loci array (all should be one because we base the sequences on one reference plus nuc conversions)
        let loci_array = Array2::<i32>::ones((candidate_variants.len(), j));

        //convert gt and lc arrays to df
        let clade_names: Vec<_> = clade_mutations_with_parents.keys().cloned().collect();
        let genotype_df =
            convert_to_dataframe(&candidate_variants, genotypes_array, &clade_names).unwrap();
        let loci_df = convert_to_dataframe(&candidate_variants, loci_array, &clade_names).unwrap();

        //then write to vcf
        write_to_vcf(&self.output, genotype_df, loci_df, self.output_bcf)?;

        //convert mutations to sequences and write sequences for each clade to fasta
        self.write_sequences_to_fasta(&clade_mutations_with_parents, &reference_path)?;
        Ok(())
    }
    fn write_sequences_to_fasta(
        &mut self,
        clade_mutations: &HashMap<String, Vec<String>>,
        reference_path: &PathBuf,
    ) -> Result<()> {
        //create a subdirectory for writing the sequences
        let outdir = &mut self.output.clone();
        outdir.push("sequences");
        dbg!(&outdir);
        fs::create_dir_all(&outdir);

        //load genome into memory
        let reference_genome = fasta::Reader::from_file(reference_path).unwrap();

        //prepare the seq
        let ref_seq_nth = reference_genome.records().nth(0).unwrap().unwrap();
        let ref_seq = ref_seq_nth.seq().to_vec();

        //iterate over candidate variants to create fasta for each clade with corresponding mutations
        // let ref_seq = reference_genome.fetch_seq_string().unwrap();
        clade_mutations.iter().for_each(|(clade, mutations)| {
            //create the file and handle
            let file = fs::File::create(format!("{}.fasta", outdir.join(clade.clone()).display()))
                .unwrap();
            let handle = io::BufWriter::new(file);
            let mut writer = bio::io::fasta::Writer::new(handle);

            //change the corresponding positions in the seq
            let mut clade_seq = ref_seq.clone();
            mutations.iter().for_each(|m| {
                let splitted: Vec<&str> = m.split('_').collect();
                let pos = splitted[0].parse::<usize>().unwrap();
                let alt_base = splitted[1].as_bytes();
                clade_seq[pos - 1] = alt_base[0];
            });

            //then write the record to file
            let record = bio::io::fasta::Record::with_attrs(clade, Some(&""), &clade_seq);
            writer
                .write_record(&record)
                .ok()
                .expect("Error writing record.");
        });
        Ok(())
    }
}

//accession_to_lineage simply converts accession ids that represent each lineage to lineage names as well as nextstrain clade
fn accession_to_lineage(
    df: &DataFrame,
    accession_to_clade_map: &HashMap<&str, &str>,
    output_path_to_mapping: &PathBuf,
) -> Result<DataFrame> {
    //clone the input df
    let mut renamed_df = df.clone();

    //rename accessions to clades and collect clade-accession mapping to a hashmap
    let mut mapping: HashMap<&str, &str> = HashMap::new();
    accession_to_clade_map.iter().for_each(|(acc, cld)| {
        let mut rename_str = format!("no clade");
        rename_str = format!("{}", cld.clone());
        renamed_df.rename(acc.clone(), &rename_str);
        mapping.insert(cld.clone(), acc.clone());
    });
    dbg!(&renamed_df.shape());

    //write mapping to file
    let mut wtr = csv::Writer::from_path(output_path_to_mapping)?;
    wtr.write_record(&vec!["Nexstrain_clade", "Accession"])?;
    mapping
        .iter()
        .for_each(|(cld, acc)| wtr.write_record(&vec![cld, acc]).unwrap());

    Ok(renamed_df)
}

fn download_resources(outdir: &PathBuf) -> Result<(PathBuf, PathBuf)> {
    let ref_link = &"https://raw.githubusercontent.com/nextstrain/ncov/1f7265a7f4e51147a38f738e3f2bcb77b9d35287/defaults/reference_seq.fasta";
    let clades_link = &"https://raw.githubusercontent.com/nextstrain/ncov/1f7265a7f4e51147a38f738e3f2bcb77b9d35287/defaults/clades.tsv";

    let ref_path = outdir.join("reference.fasta");
    let clades_path = outdir.join("clades.tsv");

    let _download_ref = {
        Command::new("wget")
            .arg("-c")
            .arg(&ref_link)
            .arg("-O")
            .arg(&ref_path)
            .status()
            .expect("failed to execute the download process for reference")
    };

    let _download_clades = {
        Command::new("wget")
            .arg("-c")
            .arg(&clades_link)
            .arg("-O")
            .arg(&clades_path)
            .status()
            .expect("failed to execute the download process for clades.tsv")
    };

    Ok((ref_path, clades_path))
}

#[derive(Debug, serde::Deserialize)]
struct Record {
    clade: String,
    gene: String,
    site: String,
    alt: Option<String>,
}

fn convert_to_dataframe(
    candidate_variants: &BTreeMap<(String, usize, String, String), Vec<usize>>,
    array: Array2<i32>,
    haplotype_names: &Vec<String>,
) -> Result<DataFrame> {
    //first initialize the dataframes with index columns that contain variant information.
    let mut genotype_df: DataFrame = df!("Index" => candidate_variants
    .iter()
    .map(|((chrom, pos, ref_base, alt_base), _)| {
        format!(
            "{},{},{},{}",
            chrom, pos, ref_base, alt_base
        )
    })
    .collect::<Series>())?;

    //construct the genotype dataframe
    for (index, haplotype_name) in haplotype_names.iter().enumerate() {
        genotype_df.with_column(Series::new(
            haplotype_name.as_str(),
            array.column(index).to_vec(),
        ))?;
    }

    //insert an ID column as many as the number of rows, right after the Index column
    genotype_df.insert_at_idx(
        1,
        Series::new("ID", Vec::from_iter(0..genotype_df.shape().0 as i64)),
    )?;
    Ok(genotype_df)
}
