use anyhow::Result;
use derive_builder::Builder;
use std::path::PathBuf;
use std::process::Command;

#[allow(dead_code)]
#[derive(Builder, Clone)]
pub struct Caller {
    alleles: PathBuf,
    genome: PathBuf,
}
impl Caller {
    pub fn call(&self) -> Result<()> {
        //TODO: generation of candidate variants from the alignment.
        Ok(())
    }
    #[allow(dead_code)]
    fn alignment(&self) -> Result<()> {
        let genome_name = format!(
            "{}{}",
            self.genome.file_name().unwrap().to_str().unwrap(),
            ".mmi"
        );
        let index = {
            Command::new("minimap2")
                .arg("-d")
                .arg(&genome_name)
                .arg(self.genome.clone())
                .status()
                .expect("failed to execute indexing process")
        };
        println!("indexing process finished with: {}", index);

        let align = {
            Command::new("minimap2")
                .args(["-a", "-t", "36"])
                .arg(&genome_name)
                .arg(self.alleles.clone())
                .status()
                .expect("failed to execute alignment process")
        };
        println!("alignment process finished with: {}", align);
        Ok(())
    }
}
