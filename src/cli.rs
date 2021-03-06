use crate::calling;
use crate::candidate_variants;
use anyhow::Result;
use rust_htslib::bcf;
use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt, Clone)]
#[structopt(
    name = "orthanq",
    about = "A haplotype caller for HLA typing and/or viral strain quantification.",
    usage = "orthanq --haplotype-counts counts.hdf5 --haplotype-variants variants.vcf \
     --haplotype-calls calls.bcf --observations observations.bcf --k-reads 20 \
     --min-norm-counts 0.01 --max-haplotypes 3 --use-evidence both --output results.tsv",
    setting = structopt::clap::AppSettings::ColoredHelp,
)]
pub enum Orthanq {
    #[structopt(
        name = "candidates",
        about = "Candidate variants are generated.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Candidates {
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
    },
    #[structopt(
        name = "call",
        about = "Call haplotypes.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Call {
        #[structopt(
            parse(from_os_str),
            long = "haplotype-counts",
            required = true,
            help = "HDF5 haplotype counts calculated by Kallisto."
        )]
        haplotype_counts: PathBuf,
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
        haplotype_calls: PathBuf,
        #[structopt(
            parse(from_os_str),
            long = "observations",
            required = true,
            help = "Variant observations by Varlociraptor."
        )]
        observations: PathBuf,
        #[structopt(
            default_value = "0.01",
            long = "min-norm-counts",
            help = "Minimum value for normalized Kallisto counts."
        )]
        min_norm_counts: f64,
        #[structopt(
            default_value = "2", //for ploidy = 2
            long = "max-haplotypes",
            help = "Expected maximum number of haplotype."
        )]
        max_haplotypes: i64,
        #[structopt(
            long,
            help = "Folder to store quality control plots for the inference of a CDF from Kallisto bootstraps for each haplotype of interest."
        )]
        output: Option<PathBuf>,
        #[structopt(long, help = "Use only kallisto evidence. (for debugging purposes)")]
        use_evidence: String,
        #[structopt(
            default_value = "20",
            long = "k-reads",
            help = "Plausible haplotypes are only those which are backed by at least k."
        )]
        k_reads: usize,
    },
}

pub fn run(opt: Orthanq) -> Result<()> {
    let opt_clone = opt.clone();
    match opt_clone {
        Orthanq::Call {
            haplotype_counts,
            haplotype_variants,
            haplotype_calls,
            observations,
            min_norm_counts,
            max_haplotypes,
            output,
            use_evidence,
            k_reads,
        } => {
            let mut caller = calling::haplotypes::CallerBuilder::default()
                .hdf5_reader(hdf5::File::open(&haplotype_counts)?)
                .haplotype_variants(bcf::Reader::from_path(&haplotype_variants)?)
                .haplotype_calls(bcf::Reader::from_path(&haplotype_calls)?)
                .min_norm_counts(min_norm_counts)
                .max_haplotypes(max_haplotypes)
                .observations(bcf::Reader::from_path(observations)?)
                .outcsv(output)
                .use_evidence(use_evidence)
                .k_reads(k_reads)
                .build()
                .unwrap();
            caller.call()?;
            Ok(())
        }
        Orthanq::Candidates { alleles, genome } => {
            let caller = candidate_variants::CallerBuilder::default()
                .alleles(alleles)
                .genome(genome)
                .build()
                .unwrap();
            caller.call()?;
            Ok(())
        }
    }
}
