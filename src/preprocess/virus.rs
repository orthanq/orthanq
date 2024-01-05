use anyhow::Result;
use derive_builder::Builder;

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;
use std::process::Stdio;

use rust_htslib::bcf;
use std::path::Path;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    genome: PathBuf,
    reads: Vec<PathBuf>,
    haplotype_variants: PathBuf,
    output: PathBuf,
}

impl Caller {
    pub fn call(&self) -> Result<()> {
        //1) bgzip and tabix the candidates vcf then, perform vg autoindex (maybe add this part to the candidates virus subcommand.)
        let thread_number = "10".to_string();

        //create the output file name in temp directory

        //get file name
        // let file_name = &self.output.file_stem().unwrap().to_str().unwrap();
        let file_name = Path::new(&self.haplotype_variants)
            .file_name()
            .unwrap()
            .to_str()
            .unwrap();
        dbg!(&file_name);

        //bgzip

        let outdir = &self.output;

        //create the folder first if it doesn't exist
        fs::create_dir_all(&outdir)?;

        let bgzip_dir = outdir.join(format!("{}.gz", file_name));
        println!("{}", bgzip_dir.display());

        let bgzip = {
            Command::new("bgzip")
                .arg("-c")
                .arg(&self.haplotype_variants)
                .stdout(Stdio::piped())
                .spawn()
                .expect("failed to execute the sorting process")
        };
        println!("Bgzip was exited with: {:?}", bgzip);

        let output = bgzip
            .wait_with_output()
            .expect("Bgzip: Failed to read stdout");

        let mut bgzipped_file = std::fs::File::create(bgzip_dir.clone())?;
        bgzipped_file.write_all(&output.stdout)?; //write with bam writer
        bgzipped_file.flush()?;

        //tabix

        let tabix = {
            Command::new("tabix")
                .arg("-p")
                .arg("vcf")
                .arg(bgzip_dir.clone())
                .stdout(Stdio::piped())
                .spawn()
                .expect("failed to execute the tabix process")
        };
        println!("Tabix was exited with: {:?}", tabix);

        let output = tabix
            .wait_with_output()
            .expect("Tabix: Failed to read stdout");

        //vg indexing

        //path to the index directory
        let giraffe_index = "resources/hprc-v1.0-mc-grch38.xg";
        //vg autoindex --workflow giraffe -r {input.genome} -v {input.variants} -p results/vg/autoindex/idx -t {threads}

        //create the output file name
        let idx_dir = outdir.join(&"idx.giraffe.gbz");
        println!("{}", bgzip_dir.display());

        let vg_index = {
            Command::new("vg")
                .arg("autoindex")
                .arg("--workflow")
                .arg("giraffe")
                .arg("-r")
                .arg(&self.genome)
                .arg("-v")
                .arg(bgzip_dir.clone())
                .arg("-p")
                .arg(idx_dir)
                .arg("-t")
                .arg(thread_number)
                .status()
                .expect("failed to execute the vg giraffe process")
        };
        println!("vg indexing was exited with: {:?}", vg_index);

        //2) vg giraffe, sorting, indexing and varlociraptor preprocess-call steps.
        Ok(())
    }
}
