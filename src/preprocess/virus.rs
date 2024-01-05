use anyhow::Result;
use derive_builder::Builder;

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;
use std::process::Stdio;

use rust_htslib::bcf;
use std::path::Path;

use tempfile::tempdir;

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

        let mut bgzipped_file = std::fs::File::create(&bgzip_dir)?;
        bgzipped_file.write_all(&output.stdout)?; //write with bam writer
        bgzipped_file.flush()?;

        //tabix

        let tabix = {
            Command::new("tabix")
                .arg("-p")
                .arg("vcf")
                .arg(&bgzip_dir)
                .stdout(Stdio::piped())
                .spawn()
                .expect("failed to execute the tabix process")
        };
        println!("Tabix was exited with: {:?}", tabix);

        let output = tabix
            .wait_with_output()
            .expect("Tabix: Failed to read stdout");

        //vg indexing

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
                .arg(&bgzip_dir)
                .arg("-p")
                .arg(&idx_dir)
                .arg("-t")
                .arg(&thread_number)
                .status()
                .expect("failed to execute the vg giraffe process")
        };
        println!("vg indexing was exited with: {:?}", vg_index);

        //2) vg giraffe, sorting, indexing and varlociraptor preprocess-call steps.

        // Create a directory inside of `std::env::temp_dir()`
        let temp_dir = tempdir()?;

        //find sample name of one of the fastq files from the read pair
        let sample_name_ext = Path::new(&self.reads[0])
            .file_name()
            .unwrap()
            .to_str()
            .unwrap();

        //get rid of underscores (PE reads contain them)
        let sample_name = sample_name_ext.split('_').collect::<Vec<&str>>()[0];

        //create the output file name in temp directory
        let file_aligned_pangenome = temp_dir.path().join(format!("{}_vg.bam", sample_name));

        //"vg giraffe -Z results/vg/autoindex/idx.giraffe.gbz -f {input.reads[0]} -f {input.reads[1]} --output-format BAM -t {threads}  > {output} 2> {log}"

        let align_pangenome = {
            Command::new("vg")
                .arg("giraffe")
                .arg("-Z")
                .arg(&idx_dir)
                .arg("-f")
                .arg(&self.reads[0])
                .arg("-f")
                .arg(&self.reads[0])
                .arg("--output-format")
                .arg("BAM")
                .arg("-t")
                .arg(&thread_number)
                .stdout(Stdio::piped())
                .spawn()
                .expect("failed to execute the vg giraffe process")
        };
        println!(
            "Alignment to pangenome was exited with: {:?}",
            align_pangenome
        );

        let output = align_pangenome
            .wait_with_output()
            .expect("Failed to read stdout");

        let mut vg_bam = std::fs::File::create(file_aligned_pangenome.clone())?;
        vg_bam.write_all(&output.stdout)?; //write with bam writer
        vg_bam.flush()?;

        //sort the resulting vg aligned file
        let file_vg_aligned_sorted = temp_dir.path().join(format!("{}_sorted.bam", sample_name));

        let vg_sort = {
            Command::new("samtools")
                .arg("sort")
                .arg(&file_aligned_pangenome)
                .arg("-o")
                .arg(&file_vg_aligned_sorted)
                .arg("-@")
                .arg(&thread_number)
                .arg("--write-index")
                .status()
                .expect("failed to execute the sorting process")
        };
        println!("The sorting was exited with: {}", vg_sort);
        println!("{}", file_vg_aligned_sorted.display());

        Ok(())
    }
}
