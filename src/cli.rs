use structopt::StructOpt;
use std::path::PathBuf;
use rust_htslib::bcf;
use crate::calling;
use anyhow::Result;

#[derive(Debug, StructOpt, Clone)]
#[structopt(
    name = "orthanq",
    about = "A haplotype caller for HLA typing and/or viral strain quantification.",
    usage = "orthanq call --haplotype-counts counts.hdf5 \
        --haplotype-variants variants.vcf --haplotype-calls calls.bcf --output results.tsv",
    setting = structopt::clap::AppSettings::ColoredHelp,
)]
pub struct Orthanq {
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
                default_value = "0.01",
                help = "Minimum value for normalized Kallisto counts."
            )]
    min_norm_counts: f64,
    #[structopt(
                long,
                help = "Folder to store quality control plots for the inference of a CDF from Kallisto bootstraps for each haplotype of interest."
            )]
    output: Option<PathBuf>,
}

pub fn run(opt: Orthanq) -> Result<()> {
    let opt_clone = opt.clone();
    // let Orthanq = Orthanq {
    //     haplotype_counts,
    //     haplotype_variants,
    //     haplotype_calls,
    //     min_norm_counts,
    //     output
    // };
    let mut caller = calling::haplotypes::CallerBuilder::default()
    .hdf5_reader(hdf5::File::open(&opt_clone.haplotype_counts)?)
    .haplotype_variants(bcf::Reader::from_path(&opt_clone.haplotype_variants)?)
    .haplotype_calls(bcf::Reader::from_path(&opt_clone.haplotype_calls)?)
    .min_norm_counts(opt_clone.min_norm_counts)
    .outcsv(opt_clone.output)
    .build()
    .unwrap();
    caller.call()?;
    
    Ok(())
}