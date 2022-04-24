use rust_htslib::bcf;

#[test]
fn check_haplotype_fractions_pure() {
    let mut output = std::path::PathBuf::new();
    output.push("test_output_pure.csv");
    let _ = &orthanq::calling::haplotypes::CallerBuilder::default()
        .hdf5_reader(hdf5::File::open("tests/Sample_HLA20547.h5").unwrap())
        .haplotype_variants(bcf::Reader::from_path("tests/hla-allele-variants_v3.vcf.gz").unwrap())
        .haplotype_calls(bcf::Reader::from_path("tests/Sample_HLA20547.bcf").unwrap())
        .min_norm_counts(0.05)
        .max_haplotypes(2)
        .outcsv(Some(output))
        .build()
        .unwrap()
        .call();

    //check if the haplotype is correct
    let mut rdr = csv::Reader::from_path("test_output_pure.csv").unwrap();
    let headers = rdr.headers().unwrap();
    assert_eq!(headers, vec!["density", "odds", "HLA:HLA20547"]);

    //check if the fraction of HLA20547 is 1.0
    let mut iter = rdr.records();
    if let Some(result) = iter.next() {
        let record = result.unwrap();
        assert_eq!(record[2], "1.00".to_string());
    }
}

#[test]
fn check_haplotype_fractions_5050() {
    let mut output = std::path::PathBuf::new();
    output.push("test_output_5050.csv");
    let _ = &orthanq::calling::haplotypes::CallerBuilder::default()
        .hdf5_reader(hdf5::File::open("tests/Sample_HLA20547_HLA24424.h5").unwrap())
        .haplotype_variants(bcf::Reader::from_path("tests/hla-allele-variants_v3.vcf.gz").unwrap())
        .haplotype_calls(bcf::Reader::from_path("tests/Sample_HLA20547_HLA24424.bcf").unwrap())
        .min_norm_counts(0.05)
        .max_haplotypes(2)
        .outcsv(Some(output))
        .build()
        .unwrap()
        .call();

    //check if the haplotypes are correct
    let mut rdr = csv::Reader::from_path("test_output_5050.csv").unwrap();

    //check if the fractions are 0.5 each
    let mut iter = rdr.records();
    if let Some(result) = iter.next() {
        let record = result.unwrap();
        assert_eq!(record[2], "0.50".to_string());
        assert_eq!(record[3], "0.50".to_string());
    }
}
