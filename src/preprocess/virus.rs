use anyhow::Result;
use derive_builder::Builder;

use std::io::Write;
use std::path::PathBuf;
use std::process::Command;
use std::process::Stdio;

use std::path::Path;

use tempfile::tempdir;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    candidates_folder: PathBuf,
    reads: Vec<PathBuf>,
    output: PathBuf,
    threads: String,
}

impl Caller {
    pub fn call(&self) -> Result<()> {
        //specify out dir
        let outdir = &self.output;

        //create the parent dir
        let mut parent = outdir.clone();
        parent.pop();
        fs::create_dir_all(parent)?;

        // Create a directory inside of `std::env::temp_dir()` for temporary files
        let temp_dir = tempdir()?;

        //1) bgzip and tabix the candidates vcf then, perform vg autoindex (maybe add this part to the candidates virus subcommand.)
        let cargo_dir = env!("CARGO_MANIFEST_DIR");
        let scenario = format!("{}/resources/scenarios/scenario.yaml", cargo_dir);

        //genome must have been downloaded in the candidate generation step:
        let ref_genome = self.candidates_folder.join("reference.fasta");

        //haplotype variantts must have been downloaded in the candidate generation step:
        let haplotype_variants = self.candidates_folder.join("candidates.vcf");
        dbg!(&haplotype_variants);
        //create the output file name in temp directory

        //get file name
        // let file_name = &self.output.file_stem().unwrap().to_str().unwrap();
        let file_name = Path::new(&haplotype_variants)
            .file_name()
            .unwrap()
            .to_str()
            .unwrap();
        dbg!(&file_name);

        //bgzip

        let bgzip_dir = temp_dir.path().join(format!("{}.gz", file_name));
        println!("{}", bgzip_dir.display());

        let bgzip = {
            Command::new("bgzip")
                .arg("-c")
                .arg(&haplotype_variants)
                .stdout(Stdio::piped())
                .spawn()
                .expect("failed to execute the zipping process")
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

        let _ = tabix
            .wait_with_output()
            .expect("Tabix: Failed to read stdout");

        //vg indexing

        //create the output file name
        let idx_name = temp_dir.path().join(&"idx");
        println!("{}", bgzip_dir.display());

        let vg_index = {
            Command::new("vg")
                .arg("autoindex")
                .arg("--workflow")
                .arg("giraffe")
                .arg("-r")
                .arg(&ref_genome)
                .arg("-v")
                .arg(&bgzip_dir)
                .arg("-p")
                .arg(&idx_name)
                .arg("-t")
                .arg(&self.threads)
                .status()
                .expect("failed to execute the vg giraffe process")
        };
        println!("vg indexing was exited with: {:?}", vg_index);

        //2) vg giraffe, sorting, indexing and varlociraptor preprocess-call steps.

        //find sample name of one of the fastq files from the read pair
        let sample_name_ext = Path::new(&self.reads[0])
            .file_name()
            .unwrap()
            .to_str()
            .unwrap();

        //get rid of underscores (PE reads contain them) (the following section needs to be rethought)
        let mut sample_name = "sample".to_string();
        if sample_name_ext.contains('_') {
            sample_name = sample_name_ext.split('_').collect::<Vec<&str>>()[0].to_string();
        } else if sample_name_ext.contains('.') {
            sample_name = sample_name_ext.split('.').collect::<Vec<&str>>()[0].to_string();
        }

        //create the output file name in temp directory
        let file_aligned_pangenome = temp_dir.path().join(format!("{}_vg.bam", sample_name));

        //"vg giraffe -Z results/vg/autoindex/idx.giraffe.gbz -f {input.reads[0]} -f {input.reads[1]} --output-format BAM -t {threads}  > {output} 2> {log}"

        let align_pangenome = {
            Command::new("vg")
                .arg("giraffe")
                .arg("-Z")
                .arg(format!("{}.giraffe.gbz", &idx_name.display()))
                .arg("-f")
                .arg(&self.reads[0])
                .arg("-f")
                .arg(&self.reads[1])
                .arg("--output-format")
                .arg("BAM")
                .arg("-t")
                .arg(&self.threads)
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
        // let file_vg_aligned_sorted = temp_dir.path().join(format!("{}_sorted.bam", sample_name));

        //print it to the parent folder of resulting bcf for debugging purposes.
        let mut parent = outdir.clone();
        parent.pop();
        let file_vg_aligned_sorted = parent.join(format!("{}_sorted.bam", sample_name));

        println!(
            "file_vg_aligned_sorteds: {}",
            file_vg_aligned_sorted.display()
        );

        let vg_sort = {
            Command::new("samtools")
                .arg("sort")
                .arg(&file_aligned_pangenome)
                .arg("-o")
                .arg(&file_vg_aligned_sorted)
                .arg("-@")
                .arg(&self.threads)
                .arg("--write-index")
                .status()
                .expect("failed to execute the sorting process")
        };
        println!("The sorting was exited with: {}", vg_sort);
        println!("{}", file_vg_aligned_sorted.display());

        //varlociraptor preprocess and call

        //preprocess
        //create the output file name
        let varlociraptor_prep_dir = temp_dir.path().join(format!("{}_obs.bcf", sample_name));
        println!(
            "varlociraptor_prep_dir: {}",
            varlociraptor_prep_dir.display()
        );

        let varlociraptor_prep = {
            Command::new("varlociraptor")
                .arg("preprocess")
                .arg("variants")
                .arg("--report-fragment-ids")
                .arg("--omit-mapq-adjustment")
                .arg("--atomic-candidate-variants")
                .arg("--candidates")
                .arg(&haplotype_variants)
                .arg(&ref_genome)
                .arg("--bam")
                .arg(&file_vg_aligned_sorted)
                .arg("--output")
                .arg(&varlociraptor_prep_dir)
                .status()
                .expect("failed to execute the varlociraptor preprocessing")
        };
        println!(
            "The varlociraptor preprocessing was exited with: {}",
            varlociraptor_prep
        );

        //call
        // "varlociraptor call variants --omit-strand-bias --omit-read-position-bias --omit-read-orientation-bias --omit-softclip-bias --omit-homopolymer-artifact-detection --omit-alt-locus-bias generic --obs sample={input.obs} " ##varlociraptor v5.3.0
        // "--scenario {input.scenario} > {output} 2> {log}"

        //scenario
        let varlociraptor_call = {
            Command::new("varlociraptor")
                .arg("call")
                .arg("variants")
                .arg("--output")
                .arg(&outdir)
                .arg("--omit-strand-bias")
                .arg("--omit-read-position-bias")
                .arg("--omit-read-orientation-bias")
                .arg("--omit-softclip-bias")
                .arg("--omit-homopolymer-artifact-detection")
                .arg("--omit-alt-locus-bias")
                .arg("generic")
                .arg("--obs")
                .arg(format!("sample={}", varlociraptor_prep_dir.display()))
                .arg("--scenario")
                .arg(&scenario)
                .status()
                .expect("failed to execute the varlociraptor calling process")
        };
        println!(
            "varlociraptor calling finished with exit status: {:?}",
            varlociraptor_call
        );

        //~fin
        Ok(())
    }
}
