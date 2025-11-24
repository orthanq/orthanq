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
            // required = true,
            help = "Input FASTQ reads belonging to the sample."
        )]
        reads: Option<Vec<PathBuf>>,
        #[structopt(
            long = "bam-input",
            // required = true,
            help = "Input BAM file (has to be aligned with BWA and annotated with read group information)."
        )]
        bam_input: Option<PathBuf>,
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
            default_value = "1",
            long = "threads",
            help = "Threads to use for tools used in preprocessing."
        )]
        threads: String,
        #[structopt(
            long,
            help = "Output BAM file used for variant calling (debugging purposes)."
        )]
        output_bam: bool,
    },
    Virus {
        #[structopt(
            long = "haplotype-variants",
            required = true,
            help = "Haplotype variants compared to a common reference."
        )]
        haplotype_variants: PathBuf,
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
            long = "output",
            help = "Output folder file to store preprocessed BAM file to be used as input in the calling step."
        )]
        output: PathBuf,
        #[structopt(
            default_value = "1",
            long = "threads",
            help = "Threads to use for tools used in preprocessing."
        )]
        threads: String,
        #[structopt(
            long,
            help = "Output BAM file used for variant calling (debugging purposes)."
        )]
        output_bam: bool,
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
        // #[structopt(
        //     long = "allele-freq",
        //     required = true,
        //     help = "allele frequencies for filterin purposes"
        // )]
        // allele_freq: PathBuf,
        #[structopt(
            long,
            help = "Folder to store candidate variants for each loci. Currently, only A, B, C and DQB1 are supported."
        )]
        output: PathBuf,
        #[structopt(
            default_value = "1",
            long = "threads",
            help = "Threads to use for minimap2 used candidate generation."
        )]
        threads: String,
        #[structopt(
            long,
            help = "Generate BCF output instead of the default VCF file format."
        )]
        output_bcf: bool,
        #[structopt(
            long,
            help = "Output BAM containing alignment of haplotypes against the reference genome."
        )]
        output_bam: bool,
    },
    Virus {
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
            default_value = "1",
            long = "threads",
            help = "Threads to use for minimap2 that is used in candidate generation."
        )]
        threads: String,
        #[structopt(
            long,
            help = "Generate BCF output instead of the default VCF file format."
        )]
        output_bcf: bool,
        #[structopt(
            long,
            help = "Output BAM containing alignment of haplotypes against the reference genome."
        )]
        output_bam: bool,
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
            help = "Folder to store results and diagnostic plots in json format."
        )]
        output: PathBuf,
        #[structopt(long, help = "Choose uniform, diploid or diploid-subclonal")]
        prior: String,
        // #[structopt(
        //     long,
        //     help = "If true, only common variants of considered haplotypes will be used in the model. Unusable at the moment."
        // )]
        // common_variants: bool,
        //Very importantly: lp_cutoff by default MUST be greater than 0.0 (higher cutoffs might lead to losing some haplotypes
        //that only very slightly differ by vaf at e.g. 0.01 fraction.
        #[structopt(
            default_value = "0.0",
            long,
            help = "Cutoff for linear program solutions."
        )]
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
            default_value = "5",
            help = "Number to constrain the number of haplotypes that LP finds. Currently more than 6 is not runtime-friendly for the Bayesian model."
        )]
        num_constraint_haplotypes: i32,
        #[structopt(long, help = "Output Datavzrd report for LP solution.")]
        output_lp_datavzrd: bool,
        #[structopt(
            long,
            help = "Sample to use in case of multisample BCFs. Sample name should match the sample name in the variant calls BCF."
        )]
        sample_name: Option<String>,
        #[structopt(
            parse(from_os_str),
            long = "limit-prediction",
            help = "Limit prediction to given input of HLA alleles e.g. A*01:01"
        )]
        limit_prediction: Option<PathBuf>,
    },
    Virus {
        #[structopt(
            long = "haplotype-variants",
            required = true,
            help = "Path to candidate variants."
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
            help = "Folder to store results and diagnostic plots in json format."
        )]
        output: PathBuf,
        #[structopt(long, help = "Choose uniform, diploid or diploid-subclonal")]
        prior: String,
        #[structopt(default_value = "0.01", help = "Cutoff for linear program solutions.")]
        lp_cutoff: f64,
        #[structopt(
            long,
            help = "Enable extension of haplotypes that are computed with linear program."
        )]
        extend_haplotypes: bool,
        #[structopt(
            long,
            default_value = "1",
            help = "Number of variant distances to extend haplotype list coming from the linear program."
        )]
        num_extend_haplotypes: i64,
        #[structopt(
            long,
            default_value = "5",
            help = "Number to constrain the number of haplotypes that LP finds."
        )]
        num_constraint_haplotypes: i32,
        #[structopt(long, help = "Output Datavzrd report for LP solution.")]
        output_lp_datavzrd: bool,
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
                output,
                prior,
                // common_variants,
                lp_cutoff,
                enable_equivalence_class_constraint,
                extend_haplotypes,
                threshold_equivalence_class,
                num_extend_haplotypes,
                num_constraint_haplotypes,
                output_lp_datavzrd,
                sample_name,
                limit_prediction,
            } => {
                let mut caller = calling::haplotypes::hla::CallerBuilder::default()
                    .haplotype_variants(bcf::Reader::from_path(haplotype_variants)?)
                    .variant_calls(bcf::Reader::from_path(variant_calls)?)
                    .xml(xml)
                    // .max_haplotypes(max_haplotypes)
                    // .min_norm_counts(min_norm_counts)
                    .output_folder(output)
                    .prior(prior)
                    // .common_variants(common_variants)
                    .lp_cutoff(lp_cutoff)
                    .enable_equivalence_class_constraint(enable_equivalence_class_constraint)
                    .extend_haplotypes(extend_haplotypes)
                    .threshold_equivalence_class(threshold_equivalence_class)
                    .num_extend_haplotypes(num_extend_haplotypes)
                    .num_constraint_haplotypes(num_constraint_haplotypes)
                    .output_lp_datavzrd(output_lp_datavzrd)
                    .sample_name(sample_name)
                    .limit_prediction(limit_prediction)
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
                extend_haplotypes,
                num_extend_haplotypes,
                num_constraint_haplotypes,
                output_lp_datavzrd,
            } => {
                let mut caller = calling::haplotypes::virus::CallerBuilder::default()
                    .haplotype_variants(bcf::Reader::from_path(haplotype_variants)?)
                    .variant_calls(bcf::Reader::from_path(variant_calls)?)
                    .output_folder(output)
                    .prior(prior)
                    .lp_cutoff(lp_cutoff)
                    .extend_haplotypes(extend_haplotypes)
                    .num_extend_haplotypes(num_extend_haplotypes)
                    .num_constraint_haplotypes(num_constraint_haplotypes)
                    .output_lp_datavzrd(output_lp_datavzrd)
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
                // allele_freq,
                output,
                threads,
                output_bcf,
                output_bam,
            } => {
                let caller = candidates::hla::CallerBuilder::default()
                    .alleles(alleles)
                    .genome(genome)
                    .xml(xml)
                    // .allele_freq(allele_freq)
                    .output(output)
                    .threads(threads)
                    .output_bcf(output_bcf)
                    .output_bam(output_bam)
                    .build()
                    .unwrap();
                caller.call()?;
                Ok(())
            }
            CandidatesKind::Virus {
                genome,
                lineages,
                output,
                threads,
                output_bcf,
                output_bam,
            } => {
                let mut caller = candidates::virus::generic::CallerBuilder::default()
                    .genome(genome)
                    .lineages(lineages)
                    .output(output)
                    .threads(threads)
                    .output_bcf(output_bcf)
                    .output_bam(output_bam)
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
                bam_input,
                haplotype_variants,
                vg_index,
                bwa_index,
                output,
                threads,
                output_bam,
            } => {
                preprocess::hla::CallerBuilder::default()
                    .genome(genome)
                    .vg_index(vg_index)
                    .bwa_index(bwa_index)
                    .reads(reads)
                    .bam_input(bam_input)
                    .haplotype_variants(haplotype_variants)
                    .output(output)
                    .threads(threads)
                    .output_bam(output_bam)
                    .build()
                    .unwrap()
                    .call()?;
                Ok(())
            }
            PreprocessKind::Virus {
                haplotype_variants,
                genome,
                reads,
                output,
                threads,
                output_bam,
            } => {
                preprocess::virus::CallerBuilder::default()
                    .haplotype_variants(haplotype_variants)
                    .genome(genome)
                    .reads(reads)
                    .output(output)
                    .threads(threads)
                    .output_bam(output_bam)
                    .build()
                    .unwrap()
                    .call()?;
                Ok(())
            }
        },
    }
}
