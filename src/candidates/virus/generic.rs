use anyhow::Result;

use derive_builder::Builder;




use crate::candidates::hla::alignment;
use crate::candidates::hla::find_variants_from_cigar;

use crate::candidates::virus::sarscov2::write_to_vcf;




















use std::io::Write;

use std::path::PathBuf;



#[allow(dead_code)]
#[derive(Builder, Clone)]
pub struct Caller {
    genome: PathBuf,
    lineages: PathBuf,
    output: PathBuf,
    threads: String,
}
impl Caller {
    pub fn call(&mut self) -> Result<()> {
        //align and sort
        alignment(
            &self.genome,
            &self.lineages,
            &self.threads,
            false,
            &self.output,
        )?;

        //find variants from cigar
        let (genotype_df, loci_df) =
            find_variants_from_cigar(&self.genome, &self.output.join("alignment_sorted.sam"))
                .unwrap();

        //write locus-wise vcf files.
        write_to_vcf(&self.output, genotype_df, loci_df)?;

        Ok(())
    }
}
