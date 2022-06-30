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
    usage = "orthanq --haplotype-variants variants.vcf \
     --haplotype-calls calls.bcf --observations observations.bcf --k-reads 20 \
     --max-haplotypes 3 --output results.tsv",
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
        #[structopt(long = "wes", help = "Specify the sample type.")]
        wes: bool,
        #[structopt(long = "wgs", help = "Specify the sample type (default).")]
        wgs: bool,
    },
    #[structopt(
        name = "call",
        about = "Call haplotypes.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Call {
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
            default_value = "2", //for ploidy = 2
            long = "max-haplotypes",
            help = "Expected maximum number of haplotype."
        )]
        max_haplotypes: usize,
        #[structopt(
            long,
            help = "Folder to store quality control plots for the inference of a CDF from Kallisto bootstraps for each haplotype of interest."
        )]
        output: Option<PathBuf>,
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
            haplotype_variants,
            haplotype_calls,
            observations,
            max_haplotypes,
            output,
            k_reads,
        } => {
            let mut caller = calling::haplotypes::CallerBuilder::default()
                .haplotype_variants(bcf::Reader::from_path(&haplotype_variants)?)
                .haplotype_calls(bcf::Reader::from_path(&haplotype_calls)?)
                .max_haplotypes(max_haplotypes)
                .observations(bcf::Reader::from_path(observations)?)
                .outcsv(output)
                .k_reads(k_reads)
                .build()
                .unwrap();
            caller.call()?;
            Ok(())
        }
        Orthanq::Candidates {
            alleles,
            genome,
            wes,
            wgs,
        } => {
            let caller = candidate_variants::CallerBuilder::default()
                .alleles(alleles)
                .genome(genome)
                .wes(wes)
                .wgs(wgs)
                .build()
                .unwrap();
            caller.call()?;
            Ok(())
        }
    }
}
