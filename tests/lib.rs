use assert_approx_eq::assert_approx_eq;
use rust_htslib::bcf;

#[test]
fn check_haplotype_fractions_pure() {
    let mut output = std::path::PathBuf::new();
    output.push("test_output_pure.csv");
    let _ = orthanq::calling::haplotypes::CallerBuilder::default()
        .haplotype_variants(bcf::Reader::from_path("tests/DQA1.vcf").unwrap())
        .variant_calls(bcf::Reader::from_path("tests/Sample_HLA00601-1.0.bcf").unwrap())
        .xml("tests/hla.xml".into())
        .common_variants(false)
        .outcsv(Some(output))
        .prior("diploid".to_string())
        .build()
        .unwrap()
        .call();

    //check if the haplotype is correct
    let mut rdr = csv::Reader::from_path("test_output_pure.csv").unwrap();
    let headers = rdr.headers().unwrap();
    assert_eq!(&headers[2], "DQA1*01:01:01");

    //check if the fraction of HLA20547 is 1.0
    let mut iter = rdr.records();
    if let Some(result) = iter.next() {
        let record = result.unwrap();
        let f: f64 = record[2].parse().unwrap();
        assert_approx_eq!(f, 1.0, 0.01);
    }
}

// #[test]
// #[ignore] //the model currently does not yield correct results for two samples each simulated from separate alleles with 50%.
// fn check_haplotype_fractions_5050() {
//     let mut output = std::path::PathBuf::new();
//     output.push("test_output_5050.csv");
//     let _ = orthanq::calling::haplotypes::CallerBuilder::default()
//         .haplotype_variants(bcf::Reader::from_path("tests/DQA1.vcf").unwrap())
//         .variant_calls(
//             bcf::Reader::from_path("tests/Sample_HLA00601-0.5_HLA21879-0.5.bcf").unwrap(),
//         )
//         .xml("tests/hla.xml".into())
//         .outcsv(Some(output))
//         .prior("diploid".to_string())
//         .common_variants(false)
//         .build()
//         .unwrap()
//         .call();

//     //check if the haplotypes are correct
//     let mut rdr = csv::Reader::from_path("test_output_5050.csv").unwrap();
//     let headers = rdr.headers().unwrap();
//     assert_eq!(
//         vec![&headers[2], &headers[3]],
//         vec!["DQA1*01:01:01", "DQA1*05:13"]
//     ); //both are from locus DQA1.

//     //check if the fractions are 0.5 each
//     let mut iter = rdr.records();
//     if let Some(result) = iter.next() {
//         let record = result.unwrap();
//         assert_eq!(&record[2], "0.50");
//         assert_eq!(&record[3], "0.50");
//     }
// }
