use crate::calling;
use crate::candidates;
use crate::preprocess;
use anyhow::Result;
use rust_htslib::bcf;
use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt, Clone)]
#[structopt(
    name = "orthanq",
    about = "A haplotype caller for HLA typing and/or viral strain quantification.",
    usage = "orthanq --haplotype-variants variants.vcf \
     --haplotype-calls calls.bcf --max-haplotypes 3 \
     --output results.tsv",
    setting = structopt::clap::AppSettings::ColoredHelp,
)]
pub enum Orthanq {
    #[structopt(
        name = "candidates",
        about = "Candidate variants are generated.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Candidates {
        #[structopt(subcommand)]
        kind: CandidatesKind,
    },
    #[structopt(
        name = "call",
        about = "Call haplotypes.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Call {
        #[structopt(subcommand)]
        kind: CallKind,
    },
    #[structopt(
        name = "preprocess",
        about = "Preprocess raw reads.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Preprocess {
        #[structopt(subcommand)]
        kind: PreprocessKind,
    },
}

#[derive(Debug, StructOpt, Clone)]
pub enum PreprocessKind {
    Hla {
        #[structopt(
            long = "genome",
            required = true,
            help = "Reference genome that is used during candidate generation."
        )]
        genome: PathBuf,
        #[structopt(long = "vg-index", required = true, help = "VG pangenome graph")]
        vg_index: PathBuf,
        #[structopt(long = "bwa-index", help = "bwa index")]
        bwa_index: Option<PathBuf>,
        #[structopt(
            long = "reads",
            required = true,
            help = "Input FASTQ reads belonging to the sample."
        )]
        reads: Vec<PathBuf>,
        #[structopt(
            parse(from_os_str),
            long = "haplotype-variants",
            required = true,
            help = "Haplotype variants compared to a common reference.", // TODO later, we will add a subcommand to generate this file with Varlociraptor as well
        )]
        haplotype_variants: PathBuf,
        #[structopt(
            long = "output",
            help = "Output BCF file to be used as input in the calling step."
        )]
        output: PathBuf,
        #[structopt(
            default_value = "2",
            help = "Threads to use for tools used in preprocessing."
        )]
        threads: String,
    },
    Virus {
        #[structopt(
            long = "candidates",
            required = true,
            help = "Folder that is used to create candidate variants."
        )]
        candidates: PathBuf,
        #[structopt(
            long = "genome",
            required = true,
            help = "Reference genome that is used during candidate generation. see 'reference.fasta' for SARS-CoV-2."
        )]
        genome: PathBuf,
        #[structopt(
            long = "reads",
            required = true,
            help = "Input FASTQ reads belonging to the sample."
        )]
        reads: Vec<PathBuf>,
        #[structopt(
            long = "output",
            help = "Output folder file to store preprocessed BAM file to be used as input in the calling step."
        )]
        output: PathBuf,
        #[structopt(
            default_value = "2",
            help = "Threads to use for tools used in preprocessing."
        )]
        threads: String,
    },
}

#[derive(Debug, StructOpt, Clone)]
pub enum CandidatesKind {
    Hla {
        #[structopt(
            long = "genome",
            required = true,
            help = "Reference genome to be used to align alleles (i.e. all existing HLA alleles) using minimap2."
        )]
        genome: PathBuf,
        #[structopt(
            long = "alleles",
            required = true,
            help = "All the alleles that exist for the gene of interest (e.g. HLA00001, HLA00002 .. for HLAs)"
        )]
        alleles: PathBuf,
        #[structopt(
            long = "xml",
            required = true,
            help = "xml file that is acquired from IMGT/HLA for the corresponding version"
        )]
        xml: PathBuf,
        #[structopt(
            long = "allele-freq",
            required = true,
            help = "allele frequencies for filterin purposes"
        )]
        allele_freq: PathBuf,
        #[structopt(
            long,
            help = "Folder to store quality control plots for the inference of a CDF from Kallisto bootstraps for each haplotype of interest."
        )]
        output: Option<PathBuf>,
        #[structopt(
            default_value = "2",
            help = "Threads to use for minimap2 used candidate generation."
        )]
        threads: String,
    },
    Virus {
        #[structopt(subcommand)]
        kind: CandidatesVirusMode,
    },
}
#[derive(Debug, StructOpt, Clone)]
pub enum CandidatesVirusMode {
    SARSCOV2 {
        #[structopt(long, help = "Folder to store candidate variants.")]
        output: PathBuf,
        #[structopt(
            default_value = "2",
            help = "Threads to use for minimap2 that is used in candidate generation."
        )]
        threads: String,
    },
    Generic {
        #[structopt(
            long = "genome",
            required = true,
            help = "Reference genome to align viral sequences using minimap2."
        )]
        genome: PathBuf,
        #[structopt(long, help = "Input fasta sequences of viral lineages.")]
        lineages: PathBuf,
        #[structopt(long, help = "Folder to store candidate variants.")]
        output: PathBuf,
        #[structopt(
            default_value = "2",
            help = "Threads to use for minimap2 that is used in candidate generation."
        )]
        threads: String,
    },
}

#[derive(Debug, StructOpt, Clone)]
pub enum CallKind {
    Hla {
        #[structopt(
            parse(from_os_str),
            long = "haplotype-variants",
            required = true,
            help = "Haplotype variants compared to a common reference.", // TODO later, we will add a subcommand to generate this file with Varlociraptor as well
        )]
        haplotype_variants: PathBuf,
        #[structopt(
            parse(from_os_str),
            long = "haplotype-calls",
            required = true,
            help = "Haplotype calls"
        )]
        variant_calls: PathBuf,
        #[structopt(
            long = "xml",
            required = true,
            help = "xml file that is acquired from IMGT/HLA for the corresponding version"
        )]
        xml: PathBuf,
        #[structopt(
            long,
            help = "Folder to store quality control plots for the inference of a CDF from Kallisto bootstraps for each haplotype of interest."
        )]
        output: PathBuf,
        #[structopt(long, help = "Choose uniform, diploid or diploid-subclonal")]
        prior: String,
        #[structopt(
            long,
            help = "If true, only common variants of considered haplotypes will be used in the model."
        )]
        common_variants: bool,
        #[structopt(default_value = "0.01", help = "Cutoff for linear program solutions.")]
        lp_cutoff: f64,
        #[structopt(
            long,
            help = "Enable equivalence based constrain during model exploration."
        )]
        enable_equivalence_class_constraint: bool,
        #[structopt(
            long,
            help = "Enable extension of haplotypes that are computed with linear program."
        )]
        extend_haplotypes: Option<bool>,
        #[structopt(
            long,
            default_value = "1",
            help = "Threshold for assigning equivalence classes."
        )]
        threshold_equivalence_class: usize,
        #[structopt(
            long,
            default_value = "3",
            help = "Number of variant distances to extend haplotype list coming from the linear program."
        )]
        num_extend_haplotypes: i64,
        #[structopt(
            long,
            default_value = "6",
            help = "Number to constrain the number of haplotypes that LP finds. Currently more than 6 is not runtime-friendly for the Bayesian model."
        )]
        num_constraint_haplotypes: i32,
    },
    Virus {
        #[structopt(
            long = "candidates-folder",
            required = true,
            help = "Folder that is used to create candidate variants."
        )]
        candidates_folder: PathBuf,
        #[structopt(
            parse(from_os_str),
            long = "haplotype-calls",
            required = true,
            help = "Haplotype calls"
        )]
        variant_calls: PathBuf,
        #[structopt(
            long,
            help = "Folder to store quality control plots for the inference of a CDF from Kallisto bootstraps for each haplotype of interest."
        )]
        output: PathBuf,
        #[structopt(long, help = "Choose uniform, diploid or diploid-subclonal")]
        prior: String,
        #[structopt(default_value = "0.01", help = "Cutoff for linear program solutions.")]
        lp_cutoff: f64,
        #[structopt(
            long,
            help = "Enable equivalence based constrain during model exploration."
        )]
        enable_equivalence_class_constraint: bool,
        #[structopt(
            long,
            help = "Enable extension of haplotypes that are computed with linear program."
        )]
        extend_haplotypes: Option<bool>,
        // #[structopt(
        //     long,
        //     default_value = "0.35",
        //     help = "Percent threshold for evaluated variants."
        // )]
        // threshold_considered_variants: f64,
        #[structopt(
            long,
            default_value = "2",
            help = "Threshold for assigning equivalence classes."
        )]
        threshold_equivalence_class: usize,
        #[structopt(
            long,
            default_value = "0",
            help = "Number of variant distances to extend haplotype list coming from the linear program."
        )]
        num_extend_haplotypes: i64, //larger than 0 is not yet supported.
        #[structopt(
            long,
            default_value = "6",
            help = "Number to constrain the number of haplotypes that LP finds. Currently more than 6 is not runtime-friendly for the Bayesian model."
        )]
        num_constraint_haplotypes: i32,
    },
}

#[derive(Debug, StructOpt, Clone)]
pub enum VirusKind {
    SARSCOV2 {
        #[structopt(
            long = "candidates-folder",
            required = true,
            help = "Folder that is used to create candidate variants."
        )]
        candidates_folder: PathBuf,
        #[structopt(
            parse(from_os_str),
            long = "haplotype-calls",
            required = true,
            help = "Haplotype calls"
        )]
        variant_calls: PathBuf,
        #[structopt(
            long,
            help = "Folder to store quality control plots for the inference of a CDF from Kallisto bootstraps for each haplotype of interest."
        )]
        output: PathBuf,
        #[structopt(long, help = "Choose uniform, diploid or diploid-subclonal")]
        prior: String,
        #[structopt(default_value = "0.01", help = "Cutoff for linear program solutions.")]
        lp_cutoff: f64,
        #[structopt(
            long,
            help = "Enable equivalence based constrain during model exploration."
        )]
        enable_equivalence_class_constraint: bool,
        // #[structopt(
        //     default_value = "0.5",
        //     help = "Percent threshold for evaluated variants."
        // )]
        // threshold_considered_variants: f64,
    },
    Generic {
        #[structopt(
            long = "candidates-folder",
            required = true,
            help = "Folder that is used to create candidate variants."
        )]
        candidates_folder: PathBuf,
        #[structopt(
            parse(from_os_str),
            long = "haplotype-calls",
            required = true,
            help = "Haplotype calls"
        )]
        variant_calls: PathBuf,
        #[structopt(long, help = "File path to store TSV table output.")]
        output: PathBuf,
        #[structopt(long, help = "Choose uniform, diploid or diploid-subclonal")]
        prior: String,
        #[structopt(default_value = "0.01", help = "Cutoff for linear program solutions.")]
        lp_cutoff: f64,
        #[structopt(
            long,
            help = "Enable equivalence based constrain during model exploration."
        )]
        enable_equivalence_class_constraint: bool,
        // #[structopt(
        //     default_value = "0.5",
        //     help = "Percent threshold for evaluated variants."
        // )]
        // threshold_considered_variants: f64,
    },
}

pub fn run(opt: Orthanq) -> Result<()> {
    let opt_clone = opt.clone();
    match opt_clone {
        Orthanq::Call { kind } => match kind {
            CallKind::Hla {
                haplotype_variants,
                variant_calls,
                xml,
                // max_haplotypes,
                // min_norm_counts,
                output,
                prior,
                common_variants,
                lp_cutoff,
                enable_equivalence_class_constraint,
                extend_haplotypes,
                threshold_equivalence_class,
                num_extend_haplotypes,
                num_constraint_haplotypes,
            } => {
                let mut caller = calling::haplotypes::hla::CallerBuilder::default()
                    .haplotype_variants(bcf::Reader::from_path(haplotype_variants)?)
                    .variant_calls(bcf::Reader::from_path(variant_calls)?)
                    .xml(xml)
                    // .max_haplotypes(max_haplotypes)
                    // .min_norm_counts(min_norm_counts)
                    .outcsv(output)
                    .prior(prior)
                    .common_variants(common_variants)
                    .lp_cutoff(lp_cutoff)
                    .enable_equivalence_class_constraint(enable_equivalence_class_constraint)
                    .extend_haplotypes(extend_haplotypes)
                    .threshold_equivalence_class(threshold_equivalence_class)
                    .num_extend_haplotypes(num_extend_haplotypes)
                    .num_constraint_haplotypes(num_constraint_haplotypes)
                    .build()
                    .unwrap();
                caller.call()?;
                Ok(())
            }
            CallKind::Virus {
                candidates_folder,
                variant_calls,
                output,
                prior,
                lp_cutoff,
                enable_equivalence_class_constraint,
                extend_haplotypes,
                threshold_equivalence_class,
                // threshold_considered_variants,
                num_extend_haplotypes,
                num_constraint_haplotypes,
            } => {
                let mut caller = calling::haplotypes::virus::CallerBuilder::default()
                    .candidates_folder(candidates_folder)
                    .variant_calls(bcf::Reader::from_path(variant_calls)?)
                    .outcsv(output)
                    .prior(prior)
                    .lp_cutoff(lp_cutoff)
                    .enable_equivalence_class_constraint(enable_equivalence_class_constraint)
                    .extend_haplotypes(extend_haplotypes)
                    .threshold_equivalence_class(threshold_equivalence_class)
                    // .threshold_considered_variants(threshold_considered_variants)
                    .num_extend_haplotypes(num_extend_haplotypes)
                    .num_constraint_haplotypes(num_constraint_haplotypes)
                    .build()
                    .unwrap();
                caller.call()?;
                Ok(())
            }
        },
        Orthanq::Candidates { kind } => match kind {
            CandidatesKind::Hla {
                alleles,
                genome,
                xml,
                allele_freq,
                output,
                threads,
            } => {
                let caller = candidates::hla::CallerBuilder::default()
                    .alleles(alleles)
                    .genome(genome)
                    .xml(xml)
                    .allele_freq(allele_freq)
                    .output(output)
                    .threads(threads)
                    .build()
                    .unwrap();
                caller.call()?;
                Ok(())
            }
            CandidatesKind::Virus { kind } => match kind {
                CandidatesVirusMode::Generic {
                    genome,
                    lineages,
                    output,
                    threads,
                } => {
                    let mut caller = candidates::virus::generic::CallerBuilder::default()
                        .genome(genome)
                        .lineages(lineages)
                        .output(output)
                        .threads(threads)
                        .build()
                        .unwrap();
                    caller.call()?;
                    Ok(())
                }
                CandidatesVirusMode::SARSCOV2 { output, threads } => {
                    let mut caller = candidates::virus::sarscov2::CallerBuilder::default()
                        .output(output)
                        .threads(threads)
                        .build()
                        .unwrap();
                    caller.call()?;
                    Ok(())
                }
            },
        },
        Orthanq::Preprocess { kind } => match kind {
            PreprocessKind::Hla {
                genome,
                reads,
                haplotype_variants,
                vg_index,
                bwa_index,
                output,
                threads,
            } => {
                preprocess::hla::CallerBuilder::default()
                    .genome(genome)
                    .vg_index(vg_index)
                    .bwa_index(bwa_index)
                    .reads(reads)
                    .haplotype_variants(haplotype_variants)
                    .output(output)
                    .threads(threads)
                    .build()
                    .unwrap()
                    .call()?;
                Ok(())
            }
            PreprocessKind::Virus {
                candidates,
                genome,
                reads,
                output,
                threads,
            } => {
                preprocess::virus::CallerBuilder::default()
                    .candidates(candidates)
                    .genome(genome)
                    .reads(reads)
                    .output(output)
                    .threads(threads)
                    .build()
                    .unwrap()
                    .call()?;
                Ok(())
            }
        },
    }
}
