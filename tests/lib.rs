use rust_htslib::bcf;

#[test]
fn check_haplotype_fractions_pure() {
    let mut output = std::path::PathBuf::new();
    output.push("test_output_pure.csv");
    let _ = orthanq::calling::haplotypes::CallerBuilder::default()
        .hdf5_reader(hdf5::File::open("tests/Sample_HLA00601-1.0_DQA1_abundance.h5").unwrap())
        .haplotype_variants(bcf::Reader::from_path("tests/hla-allele-variants_v4.vcf.gz").unwrap())
        .haplotype_calls(bcf::Reader::from_path("tests/Sample_HLA00601-1.0.bcf").unwrap())
        .observations(bcf::Reader::from_path("tests/Sample_HLA00601-1.0_obs.bcf").unwrap())
        .min_norm_counts(0.01)
        .max_haplotypes(2)
        .use_evidence("both".to_string())
        .outcsv(Some(output))
        .build()
        .unwrap()
        .call();

    //check if the haplotype is correct
    let mut rdr = csv::Reader::from_path("test_output_pure.csv").unwrap();
    let headers = rdr.headers().unwrap();
    assert_eq!(&headers[2], "HLA:HLA00601");

    //check if the fraction of HLA20547 is 1.0
    let mut iter = rdr.records();
    if let Some(result) = iter.next() {
        let record = result.unwrap();
        assert_eq!(&record[2], "1.00");
    }
}

#[test]
fn check_haplotype_fractions_5050() {
    let mut output = std::path::PathBuf::new();
    output.push("test_output_5050.csv");
    let _ = orthanq::calling::haplotypes::CallerBuilder::default()
        .hdf5_reader(
            hdf5::File::open("tests/Sample_HLA00601-0.5_HLA21879-0.5_DQA1_abundance.h5").unwrap(),
        )
        .haplotype_variants(bcf::Reader::from_path("tests/hla-allele-variants_v4.vcf.gz").unwrap())
        .haplotype_calls(
            bcf::Reader::from_path("tests/Sample_HLA00601-0.5_HLA21879-0.5.bcf").unwrap(),
        )
        .observations(
            bcf::Reader::from_path("tests/Sample_HLA00601-0.5_HLA21879-0.5_obs.bcf").unwrap(),
        )
        .min_norm_counts(0.01)
        .max_haplotypes(2)
        .use_evidence("both".to_string())
        .outcsv(Some(output))
        .build()
        .unwrap()
        .call();

    //check if the haplotypes are correct
    let mut rdr = csv::Reader::from_path("test_output_5050.csv").unwrap();
    let headers = rdr.headers().unwrap();
    assert_eq!(
        vec![&headers[2], &headers[3]],
        vec!["HLA:HLA00601", "HLA:HLA21879"]
    ); //both are from locus DQA1.

    //check if the fractions are 0.5 each
    let mut iter = rdr.records();
    if let Some(result) = iter.next() {
        let record = result.unwrap();
        assert_eq!(&record[2], "0.50");
        assert_eq!(&record[3], "0.50");
    }
}
