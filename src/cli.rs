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
     --haplotype-calls calls.bcf --max-haplotypes 3 \
     --use-evidence both --output results.tsv",
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
        #[structopt(
            long = "xml",
            required = true,
            help = "xml file that is acquired from IMGT/HLA for the corresponding version"
        )]
        xml: PathBuf,
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
        variant_calls: PathBuf,
        // #[structopt(
        //     parse(from_os_str),
        //     long = "observations",
        //     required = true,
        //     help = "Variant observations by Varlociraptor."
        // )]
        // observations: PathBuf,
        #[structopt(
            default_value = "5",
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
        #[structopt(long, help = "Choose uniform or diploid")]
        prior: String,
    },
}

pub fn run(opt: Orthanq) -> Result<()> {
    let opt_clone = opt.clone();
    match opt_clone {
        Orthanq::Call {
            haplotype_variants,
            variant_calls,
            max_haplotypes,
            output,
            use_evidence,
            prior
        } => {
            let mut caller = calling::haplotypes::CallerBuilder::default()
                .haplotype_variants(bcf::Reader::from_path(&haplotype_variants)?)
                .variant_calls(bcf::Reader::from_path(&variant_calls)?)
                .max_haplotypes(max_haplotypes)
                .outcsv(output)
                .use_evidence(use_evidence)
                .prior(prior)
                .build()
                .unwrap();
            caller.call()?;
            Ok(())
        }
        Orthanq::Candidates {
            alleles,
            genome,
            xml,
            wes,
            wgs,
            output,
        } => {
            let caller = candidate_variants::CallerBuilder::default()
                .alleles(alleles)
                .genome(genome)
                .xml(xml)
                .wes(wes)
                .wgs(wgs)
                .output(output)
                .build()
                .unwrap();
            caller.call()?;
            Ok(())
        }
    }
}
