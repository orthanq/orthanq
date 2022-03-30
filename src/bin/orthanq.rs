use anyhow::Result;
use orthanq::cli::{run, Orthanq};
use structopt::StructOpt;

pub(crate) fn main() -> Result<()> {
    let opt = Orthanq::from_args();

    // setup logger
    fern::Dispatch::new()
        .level(log::LevelFilter::Info)
        .chain(std::io::stderr())
        .apply()
        .unwrap();

    run(opt)?;

    Ok(())
}
