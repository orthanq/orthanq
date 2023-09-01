use rust_htslib::bcf;

#[test]
fn check_haplotype_fractions_5050() {
    let mut output = std::path::PathBuf::new();
    output.push("test_output.csv");
    let _ = orthanq::calling::haplotypes::CallerBuilder::default()
        .haplotype_variants(bcf::Reader::from_path("tests/B.vcf").unwrap())
        .variant_calls(
            bcf::Reader::from_path("tests/Sample_HLA00318-0.5_HLA00319-0.5_B.bcf").unwrap(),
        )
        .xml("tests/hla.xml".into())
        .common_variants(false)
        .outcsv(Some(output))
        .prior("diploid".to_string())
        .build()
        .unwrap()
        .call();

    //check if the haplotype is correct
    let mut rdr = csv::Reader::from_path("test_output.csv").unwrap();

    //check if the fractions are 0.5 each
    let mut iter = rdr.records();
    if let Some(result) = iter.next() {
        let record = result.unwrap();
        assert_eq!(&record[2], "0.50"); //44:02:01
        assert_eq!(&record[8], "0.50"); //44:03:01
    }
}
