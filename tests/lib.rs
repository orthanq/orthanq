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
        .lp_cutoff(0.01)
        .build()
        .unwrap()
        .call();

    //check if the haplotype is correct
    let mut rdr = csv::Reader::from_path("test_output.csv").unwrap();

    //access haplotype names
    let headers = rdr.headers().unwrap().clone();

    //check if the fractions are 0.5 each
    let mut iter = rdr.records();
    if let Some(result) = iter.next() {
        let record = result.unwrap();
        let mut check_two_alleles = 0;
        for (fraction,header) in record.iter().zip(headers.iter()){
            if header.contains("B*"){
                let splitted = header.split(':').collect::<Vec<&str>>();
                let first_two = splitted[0].to_owned() + ":" + splitted[1];
                if (first_two == "B*44:02" && fraction == "0.50") || (first_two == "B*44:03" && fraction == "0.50"){
                    check_two_alleles += 1;
                }   
            }
        }
        let mut check = false;
        if check_two_alleles == 2 {
            check = true
        }
        assert!(check);
    }
}
