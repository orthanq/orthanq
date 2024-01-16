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
            long = "output BCF",
            help = "Output BCF file to be used as input in the calling step."
        )]
        output: PathBuf,
    },
    Virus {
        #[structopt(
            long = "genome",
            required = true,
            help = "Reference genome that is used during candidate generation."
        )]
        genome: PathBuf,
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
            help = "Output folder file to store preprocessed BAM file to be used as input in the calling step."
        )]
        output: PathBuf,
    },
}

#[derive(Debug, StructOpt, Clone)]
pub enum CandidatesKind {
    Hla {
        #[structopt(
            long = "genome",
            required = true,
            help = "Reference genome to be used to align alleles or viral sequences (e.g. all existing HLA alleles) using minimap2."
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
        #[structopt(long = "wes", help = "Specify the sample type.")]
        wes: bool,
        #[structopt(long = "wgs", help = "Specify the sample type (default).")]
        wgs: bool,
        #[structopt(
            long,
            help = "Folder to store quality control plots for the inference of a CDF from Kallisto bootstraps for each haplotype of interest."
        )]
        output: Option<PathBuf>,
    },
    Virus {
        #[structopt(
            long = "genome",
            required = true,
            help = "Reference genome to be used to align alleles or viral sequences (e.g. all existing HLA alleles) using minimap2."
        )]
        genome: PathBuf,
        #[structopt(
            long = "alleles",
            required = true,
            help = "All the alleles that exist for the gene of interest (e.g. HLA00001, HLA00002 .. for HLAs)"
        )]
        alleles: PathBuf,
        #[structopt(
            long,
            help = "Folder to store quality control plots for the inference of a CDF from Kallisto bootstraps for each haplotype of interest."
        )]
        output: Option<PathBuf>,
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
        #[structopt(long, help = "Enable equivalence based constrain during model exploration.")]
        enable_equivalence_class_constraint: bool
    },
    Virus {
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
            long,
            help = "Folder to store quality control plots for the inference of a CDF from Kallisto bootstraps for each haplotype of interest."
        )]
        output: PathBuf,
        #[structopt(long, help = "Choose uniform, diploid or diploid-subclonal")]
        prior: String,
        #[structopt(default_value = "0.01", help = "Cutoff for linear program solutions.")]
        lp_cutoff: f64,
        #[structopt(long, help = "Enable equivalence based constrain during model exploration.")]
        enable_equivalence_class_constraint: bool
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
                enable_equivalence_class_constraint
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
                    .build()
                    .unwrap();
                caller.call()?;
                Ok(())
            }
            CallKind::Virus {
                haplotype_variants,
                variant_calls,
                output,
                prior,
                lp_cutoff,
                enable_equivalence_class_constraint
            } => {
                let mut caller = calling::haplotypes::virus::CallerBuilder::default()
                    .haplotype_variants(bcf::Reader::from_path(haplotype_variants)?)
                    .variant_calls(bcf::Reader::from_path(variant_calls)?)
                    .outcsv(output)
                    .prior(prior)
                    .lp_cutoff(lp_cutoff)
                    .enable_equivalence_class_constraint(enable_equivalence_class_constraint)
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
                wes,
                wgs,
                output,
            } => {
                let caller = candidates::hla::CallerBuilder::default()
                    .alleles(alleles)
                    .genome(genome)
                    .xml(xml)
                    .allele_freq(allele_freq)
                    .wes(wes)
                    .wgs(wgs)
                    .output(output)
                    .build()
                    .unwrap();
                caller.call()?;
                Ok(())
            }
            CandidatesKind::Virus {
                alleles,
                genome,
                output,
            } => {
                let caller = candidates::virus::CallerBuilder::default()
                    .alleles(alleles)
                    .genome(genome)
                    .output(output)
                    .build()
                    .unwrap();
                caller.call()?;
                Ok(())
            }
        },
        Orthanq::Preprocess { kind } => match kind {
            PreprocessKind::Hla {
                genome,
                reads,
                haplotype_variants,
                output,
            } => {
                preprocess::hla::CallerBuilder::default()
                    .genome(genome)
                    .reads(reads)
                    .haplotype_variants(haplotype_variants)
                    .output(output)
                    .build()
                    .unwrap()
                    .call()?;
                Ok(())
            }
            PreprocessKind::Virus {
                genome,
                reads,
                haplotype_variants,
                output,
            } => {
                preprocess::virus::CallerBuilder::default()
                    .genome(genome)
                    .reads(reads)
                    .haplotype_variants(haplotype_variants)
                    .output(output)
                    .build()
                    .unwrap()
                    .call()?;
                Ok(())
            }
        },
    }
}
