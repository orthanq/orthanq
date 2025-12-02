# Changelog

## [1.20.0](https://github.com/orthanq/orthanq/compare/v1.19.0...v1.20.0) (2025-12-02)


### Features

* enforce candidate vcf to given list of input alleles for debugging puposes ([#103](https://github.com/orthanq/orthanq/issues/103)) ([d5e3c8a](https://github.com/orthanq/orthanq/commit/d5e3c8a32a3146af452f27854a56f48e50d38dda))
* output 4-field resolution ([#102](https://github.com/orthanq/orthanq/issues/102)) ([fc7ea28](https://github.com/orthanq/orthanq/commit/fc7ea28ffc22de2c87b87be25a45566c180acb52))

## [1.19.0](https://github.com/orthanq/orthanq/compare/v1.18.1...v1.19.0) (2025-11-24)


### Features

* add support for alleles DQA1 and DRB1 for candidate variant generation ([#109](https://github.com/orthanq/orthanq/issues/109)) ([690c6b1](https://github.com/orthanq/orthanq/commit/690c6b1a07a83778b1fbf4e652b8c66fed32ef9b))

## [1.18.1](https://github.com/orthanq/orthanq/compare/v1.18.0...v1.18.1) (2025-11-19)


### Bug Fixes

* add coordindates of other loci ([#104](https://github.com/orthanq/orthanq/issues/104)) ([e51f387](https://github.com/orthanq/orthanq/commit/e51f3872066bab83871538bf86f10c071ba434a2))
* do not create datavzrd config yaml at all times ([#106](https://github.com/orthanq/orthanq/issues/106)) ([3317a54](https://github.com/orthanq/orthanq/commit/3317a5477fe76ca33280a6325e2722d4ca768337))
* do not output eq graph by default ([#107](https://github.com/orthanq/orthanq/issues/107)) ([02c4c27](https://github.com/orthanq/orthanq/commit/02c4c27e52f12956ba0f3fb31f80f862b4e7c1d8))

## [1.18.0](https://github.com/orthanq/orthanq/compare/v1.17.0...v1.18.0) (2025-10-23)


### Features

* arrow plot ([#79](https://github.com/orthanq/orthanq/issues/79)) ([d1c65cc](https://github.com/orthanq/orthanq/commit/d1c65cc9d4b55d036d05d03b7867e063f0d785b3))


### Bug Fixes

* add prob of variant to arrow plot through opacity ([#99](https://github.com/orthanq/orthanq/issues/99)) ([a36d2db](https://github.com/orthanq/orthanq/commit/a36d2db0e6c191b80b380811601ada173fb41917))
* exclude variants that arent on the locus ([#101](https://github.com/orthanq/orthanq/issues/101)) ([695e805](https://github.com/orthanq/orthanq/commit/695e8057c549924b30773734beaa7fe313820be8))
* output only nonzero fraction haplotypes & display variants that are present in one of the resulting haplotypes ([#97](https://github.com/orthanq/orthanq/issues/97)) ([56f9383](https://github.com/orthanq/orthanq/commit/56f938320ce971b887644f28727b9c0e3cee4bee))
* sort by chromosomal locations ([#98](https://github.com/orthanq/orthanq/issues/98)) ([caa0c9f](https://github.com/orthanq/orthanq/commit/caa0c9f6a03ae95da914a243aac65835560bbfa0))
* use filled circles when the expectation meets the predicted fractions ([#100](https://github.com/orthanq/orthanq/issues/100)) ([f11b607](https://github.com/orthanq/orthanq/commit/f11b6076d808bc9e5428a5350e393ae1f4fdc69e))

## [1.17.0](https://github.com/orthanq/orthanq/compare/v1.16.1...v1.17.0) (2025-10-10)


### Features

* support output to directory instead of CSV file ([#92](https://github.com/orthanq/orthanq/issues/92)) ([42bdc7f](https://github.com/orthanq/orthanq/commit/42bdc7f662f48362c09eb364bb1b26c03cf2b440))

### [1.16.1](https://www.github.com/orthanq/orthanq/compare/v1.16.0...v1.16.1) (2025-10-07)


### Bug Fixes

* Update haplotype bar plot config with view and axis changes ([#88](https://www.github.com/orthanq/orthanq/issues/88)) ([b613f3e](https://www.github.com/orthanq/orthanq/commit/b613f3e7732f320b83354073b4d073a784abd7a2))

## [1.16.0](https://www.github.com/orthanq/orthanq/compare/v1.15.0...v1.16.0) (2025-10-01)


### Features

* --sample-name parameter for hla calling ([#86](https://www.github.com/orthanq/orthanq/issues/86)) ([e17a5c4](https://www.github.com/orthanq/orthanq/commit/e17a5c4445608e8e2b11ad300dfc2c624897cb3d))


### Bug Fixes

* handle missing DP ([#85](https://www.github.com/orthanq/orthanq/issues/85)) ([7a090e2](https://www.github.com/orthanq/orthanq/commit/7a090e2815a2c2f147ab093af261e9065e900343))

## [1.15.0](https://www.github.com/orthanq/orthanq/compare/v1.14.0...v1.15.0) (2025-09-21)


### Features

* implement bcf export option for candidate variant generation ([#82](https://www.github.com/orthanq/orthanq/issues/82)) ([6f82e4b](https://www.github.com/orthanq/orthanq/commit/6f82e4bdceaa56496ac259f7367959f0eae3aba9))
* implement optional parameter for outputting bam file produced with alleles to reference genome alignment and configure virus output to the given filepath instead of the previous candidates.vcf ([#84](https://www.github.com/orthanq/orthanq/issues/84)) ([9a827ba](https://www.github.com/orthanq/orthanq/commit/9a827baa04a29fee1c1b1c72526c1532a5e3c676))

## [1.14.0](https://www.github.com/orthanq/orthanq/compare/v1.13.0...v1.14.0) (2025-07-22)


### Features

* separately align HLA alleles to genomic regions ([#80](https://www.github.com/orthanq/orthanq/issues/80)) ([1b321f6](https://www.github.com/orthanq/orthanq/commit/1b321f6e0f49b95c75f9c87c59dc1f5a0d5450aa))


### Bug Fixes

* update rust bio version ([#77](https://www.github.com/orthanq/orthanq/issues/77)) ([af7c867](https://www.github.com/orthanq/orthanq/commit/af7c867b57e8671560822ce8f0c0940a0a6267dc))

## [1.13.0](https://www.github.com/orthanq/orthanq/compare/v1.12.1...v1.13.0) (2025-05-28)


### Features

* deactivate allele frequency filtering ([#75](https://www.github.com/orthanq/orthanq/issues/75)) ([76e0409](https://www.github.com/orthanq/orthanq/commit/76e040920a9bcd160da0dde573421360b4de6497))


### Bug Fixes

* fix regex for datavzrd report ([4e44106](https://www.github.com/orthanq/orthanq/commit/4e44106f9b0ed58677fd16d8c9a772084ddf0f72))
* include size channel for tick plot and normalize values by max ([081bfee](https://www.github.com/orthanq/orthanq/commit/081bfee2b1c0e9579449fa3dd80c209b32282fff))

### [1.12.1](https://www.github.com/orthanq/orthanq/compare/v1.12.0...v1.12.1) (2025-05-12)


### Bug Fixes

* update datavzrd ([#72](https://www.github.com/orthanq/orthanq/issues/72)) ([8504aac](https://www.github.com/orthanq/orthanq/commit/8504aac31f791dbc075ef12f40d38a237492af09))

## [1.12.0](https://www.github.com/orthanq/orthanq/compare/v1.10.0...v1.12.0) (2025-05-08)


### Features

* apply the max of PROB_PRESENT or PROB_ABSENT for weighting the LP constraints
* add pinned feature to template for first 4 columns of lp datavzrd report ([#66](https://www.github.com/orthanq/orthanq/issues/66)) ([e8ee175](https://www.github.com/orthanq/orthanq/commit/e8ee175838288f3e896a9c4b7c511f04f3150357))
* export searchable and filterable datavzrd report for lp solutions ([#59](https://www.github.com/orthanq/orthanq/issues/59)) ([571ed2f](https://www.github.com/orthanq/orthanq/commit/571ed2f1ed8a90a7ee67107934b76d8c729fbb01))
* replace bubble plot with tick marks and use normal probs instead of logprobs in the final solution plot ([#65](https://www.github.com/orthanq/orthanq/issues/65)) ([dad9322](https://www.github.com/orthanq/orthanq/commit/dad93229d939057bf7166d7a5bc7c92091990cc3))
* use only probable variants for LP ([#64](https://www.github.com/orthanq/orthanq/issues/64)) ([7d35bf1](https://www.github.com/orthanq/orthanq/commit/7d35bf1cc43f94e93f10e5618f078d8813ee6998))


### Bug Fixes

* add chr name to plots ([#56](https://www.github.com/orthanq/orthanq/issues/56)) ([606f9fe](https://www.github.com/orthanq/orthanq/commit/606f9fe6f9ed135f84cc7d919571b42e4a75f641))
* find similar haplotypes by identical variant set of lp haplotypes and extend lp haplotypes accordingly and refactor VariantCalls struct ([#67](https://www.github.com/orthanq/orthanq/issues/67)) ([2acdbc0](https://www.github.com/orthanq/orthanq/commit/2acdbc0d327919cba35fe7147abb49fee800dce3))
* fix g group outputting ([#61](https://www.github.com/orthanq/orthanq/issues/61)) ([bb6038f](https://www.github.com/orthanq/orthanq/commit/bb6038f9ed53d71ac3e3247993c54ed25e437232))
* fix retrieval of confirmed field from xml file ([#58](https://www.github.com/orthanq/orthanq/issues/58)) ([4c4ed31](https://www.github.com/orthanq/orthanq/commit/4c4ed316519be5fc88b96b29ffa5cd5d42a954bb))
* output empty output when there are no calls for LP ([#69](https://www.github.com/orthanq/orthanq/issues/69)) ([ee92b9a](https://www.github.com/orthanq/orthanq/commit/ee92b9abd6abc631b4905f93bbe1ee0de5f40a2c))
* provide 2-field csv and G rgoups in empty output ([#70](https://www.github.com/orthanq/orthanq/issues/70)) ([8fffcd8](https://www.github.com/orthanq/orthanq/commit/8fffcd80680c78614834f8559c5291b6c57f8a4e))
* remove unnecessary G from the outputted name and handle case where there is no g group field for the allele in the xml file and handle haplotypes that do not have direct match with the names in the xml ([#63](https://www.github.com/orthanq/orthanq/issues/63)) ([ecc19ad](https://www.github.com/orthanq/orthanq/commit/ecc19ad3a6d1091a9087e8c95a3be82b06d3e4a1))
* use the length of haplotypes in case the length is smaller than num_constraint_haplotypes and decrement this number if the solution is infeasible ([71e6170](https://www.github.com/orthanq/orthanq/commit/71e6170f4ebdd9f2de19800b6fb395aa27a21ea7))


### Miscellaneous Chores

* release 1.12.0 ([530ba7e](https://www.github.com/orthanq/orthanq/commit/530ba7e020f28f3c2a427c3c2a25334d81d1e130))

## [1.10.0](https://www.github.com/orthanq/orthanq/compare/v1.9.3...v1.10.0) (2025-03-24)


### Features

* display nucleotide changes in plots ([#53](https://www.github.com/orthanq/orthanq/issues/53)) ([3f870fc](https://www.github.com/orthanq/orthanq/commit/3f870fc5832bd8d0964ba77d7668d2b18b6be420))

### [1.9.3](https://www.github.com/orthanq/orthanq/compare/v1.9.2...v1.9.3) (2025-03-17)


### Bug Fixes

* preprocessing: implement parameter for outputting final bam file to parent folder of the output ([#51](https://www.github.com/orthanq/orthanq/issues/51)) ([30b1be9](https://www.github.com/orthanq/orthanq/commit/30b1be9af64822aaa59ee8896c0bd8556dfb02a8))

### [1.9.2](https://www.github.com/orthanq/orthanq/compare/v1.9.1...v1.9.2) (2025-03-17)


### Bug Fixes

* fix subtraction in log space ([#47](https://www.github.com/orthanq/orthanq/issues/47)) ([9a85787](https://www.github.com/orthanq/orthanq/commit/9a85787ea81761b5280ee803a409fb94d78b284d))
* handle each process separately if they are not successful ([#48](https://www.github.com/orthanq/orthanq/issues/48)) ([787ad94](https://www.github.com/orthanq/orthanq/commit/787ad94b0bea835aea7f3306900e8f3228fb5d73))

### [1.9.1](https://www.github.com/orthanq/orthanq/compare/v1.9.0...v1.9.1) (2025-03-05)


### Bug Fixes

* fix threads usage and refactor ([#45](https://www.github.com/orthanq/orthanq/issues/45)) ([8aa406a](https://www.github.com/orthanq/orthanq/commit/8aa406af55f8384da9535d85d978e101f51e7626))

## [1.9.0](https://www.github.com/orthanq/orthanq/compare/v1.8.0...v1.9.0) (2025-03-04)


### Features

* constrain lp haplotypes to 5 by default, disable the threshold for number of variants & include all variants (nonzero DP), do not compute HD, add all variants in final solution plot ([#42](https://www.github.com/orthanq/orthanq/issues/42)) ([9926086](https://www.github.com/orthanq/orthanq/commit/9926086b61efca185e2b9840f8e659403e7690b8))

## [1.8.0](https://www.github.com/orthanq/orthanq/compare/v1.7.9...v1.8.0) (2025-02-26)


### Features

* use bam input for preprocessing hla if provided ([#43](https://www.github.com/orthanq/orthanq/issues/43)) ([0e6e1cf](https://www.github.com/orthanq/orthanq/commit/0e6e1cf91dccd7cbca6d99e8bfddbb58c547fa91))
* use lp selected haplotypes in model and extend resulting table by 0 distance haplotypes ([#39](https://www.github.com/orthanq/orthanq/issues/39)) ([b7f0561](https://www.github.com/orthanq/orthanq/commit/b7f0561fac85f48982b527841694f2353a38aff0))

### [1.7.9](https://www.github.com/orthanq/orthanq/compare/v1.7.8...v1.7.9) (2024-10-07)


### Bug Fixes

* remove if-else block from the linear interpolation func ([#35](https://www.github.com/orthanq/orthanq/issues/35)) ([d0b4562](https://www.github.com/orthanq/orthanq/commit/d0b4562049ecc1458b0e75f5914d239572610a70))

### [1.7.8](https://www.github.com/orthanq/orthanq/compare/v1.7.7...v1.7.8) (2024-08-24)


### Bug Fixes

* handle allele name parsing from xml with and without HLA prefix ([#33](https://www.github.com/orthanq/orthanq/issues/33)) ([9bcf1f9](https://www.github.com/orthanq/orthanq/commit/9bcf1f9d271ff6cbc151b51163081d30b31b434a))

### [1.7.7](https://www.github.com/orthanq/orthanq/compare/v1.7.6...v1.7.7) (2024-08-23)


### Bug Fixes

* fix threshold for considered variants to long argument and lower threshold ([#31](https://www.github.com/orthanq/orthanq/issues/31)) ([6062f8a](https://www.github.com/orthanq/orthanq/commit/6062f8a58b7e9ccaab1e1c9ff76c9a36147e6eed))

### [1.7.6](https://www.github.com/orthanq/orthanq/compare/v1.7.5...v1.7.6) (2024-08-22)


### Bug Fixes

* fix equivalence class generation for hla and virus case separately; further, implement a flag for extension of haplotype list after linear program computation ([#29](https://www.github.com/orthanq/orthanq/issues/29)) ([cc7a92d](https://www.github.com/orthanq/orthanq/commit/cc7a92df98e6f04cb63d5bf6ac68914644792772))

### [1.7.5](https://www.github.com/orthanq/orthanq/compare/v1.7.4...v1.7.5) (2024-08-19)


### Bug Fixes

* implement user defined parameters threshold for equivalence classes and extending haplotype list; also properly adjust path of equivalence class graph ([#27](https://www.github.com/orthanq/orthanq/issues/27)) ([b9f8255](https://www.github.com/orthanq/orthanq/commit/b9f8255818360b16c76abaa41bc20e326a01f2c6))

### [1.7.4](https://www.github.com/orthanq/orthanq/compare/v1.7.3...v1.7.4) (2024-08-15)


### Bug Fixes

* fixate good_lp version and do not use wildcard constraint in cargo.toml ([68afd30](https://www.github.com/orthanq/orthanq/commit/68afd305236f85893fe68cec135c79d2ef19e983))
* fix scenario path for hla preprocessing

### [1.7.3](https://www.github.com/orthanq/orthanq/compare/v1.7.2...v1.7.3) (2024-08-14)


### Bug Fixes

* add description and license info to cargo.toml ([#23](https://www.github.com/orthanq/orthanq/issues/23)) ([a3036e9](https://www.github.com/orthanq/orthanq/commit/a3036e9b265c164fccc3ef1abd672de3feefe2e2))
* fix no such file or directory error by reading scenario to string and writing it to a temp.

### [1.7.2](https://www.github.com/orthanq/orthanq/compare/v1.7.1...v1.7.2) (2024-08-14)


### Bug Fixes

* perform interpolation of allele frequency distribution on the probability scale instead of PHRED scale; further, this changes AF distribution display to a bubble chart ([#21](https://www.github.com/orthanq/orthanq/issues/21)) ([309d33a](https://www.github.com/orthanq/orthanq/commit/309d33a9f2517b9e6362e693b666dedb667e4407))

### [1.7.1](https://www.github.com/orthanq/orthanq/compare/v1.7.0...v1.7.1) (2024-08-05)


### Bug Fixes

* fix preprocess hla outputs ([#19](https://www.github.com/orthanq/orthanq/issues/19)) ([34c78ae](https://www.github.com/orthanq/orthanq/commit/34c78aebfcdd7b77b2eb9d1a051ee45cb51e2dcd))

## [1.7.0](https://www.github.com/orthanq/orthanq/compare/v1.6.0...v1.7.0) (2024-07-17)


### Features

* implement generic virus application ([#17](https://www.github.com/orthanq/orthanq/issues/17)) ([05bfe4e](https://www.github.com/orthanq/orthanq/commit/05bfe4e1f1efc90247da2db1d94e9e5b60326925))

## [1.6.0](https://www.github.com/orthanq/orthanq/compare/v1.5.0...v1.6.0) (2024-07-08)


### Features

* set a configurable threshold that data has sufficient observation ([#15](https://www.github.com/orthanq/orthanq/issues/15)) ([09bea42](https://www.github.com/orthanq/orthanq/commit/09bea426a813e1231380a030115710b318f7a048))

## [1.5.0](https://www.github.com/orthanq/orthanq/compare/v1.4.0...v1.5.0) (2024-06-25)


### Features

* optimize model by pruning search space ([#8](https://www.github.com/orthanq/orthanq/issues/8)) ([d16cbe8](https://www.github.com/orthanq/orthanq/commit/d16cbe853f61e1ec9e517632a7453aae42f5ab97))

## [1.4.0](https://www.github.com/orthanq/orthanq/compare/v1.3.0...v1.4.0) (2024-06-20)


### Features

* download nextstrain resources and generate candidate variants from nuc conversions ([#12](https://www.github.com/orthanq/orthanq/issues/12)) ([f1ffda9](https://www.github.com/orthanq/orthanq/commit/f1ffda92f328d4078746b0b374b9e9dc93f8037e))

## [1.3.0](https://www.github.com/orthanq/orthanq/compare/v1.2.0...v1.3.0) (2024-04-12)


### Features

* breaking change: orthanq virus candidates: use the oldest submitted genome to represent each nextstrain in candidate haplotypes avoiding  irrelevant pango lineages ([96f397b](https://www.github.com/orthanq/orthanq/commit/96f397b0ba99f0f695b93b75303da3ed0d269b52))

## [1.2.0](https://www.github.com/orthanq/orthanq/compare/v1.1.0...v1.2.0) (2024-02-05)


### Features

* retrieve metadata and download data packages belonging to sarscov2 and use in candidate variant generation ([#6](https://www.github.com/orthanq/orthanq/issues/6)) ([dfce2a5](https://www.github.com/orthanq/orthanq/commit/dfce2a57d257e075765230b16cbaac1cd6d00ae7))

## [1.1.0](https://www.github.com/orthanq/orthanq/compare/v1.0.1...v1.1.0) (2024-01-05)


### Features

* enable viral quantification  ([#4](https://www.github.com/orthanq/orthanq/issues/4)) ([b945239](https://www.github.com/orthanq/orthanq/commit/b945239b81293f6520ad7cb0def6cbb235dca517))

### [1.0.1](https://www.github.com/orthanq/orthanq/compare/v1.0.0...v1.0.1) (2023-12-21)


### Miscellaneous Chores

* release 1.0.1 ([af15015](https://www.github.com/orthanq/orthanq/commit/af15015df8b83c89c74ab571d2dd61582355c980))

## 1.0.0 (2023-12-21)


### Features

* output 2-field info by summing up densitites of events that have identical first two fields ([dbd192c](https://www.github.com/orthanq/orthanq/commit/dbd192c3349deda41938d346028b3cfd617f8d40))


### Miscellaneous Chores

* release 1.0.0 ([827179b](https://www.github.com/orthanq/orthanq/commit/827179b814ce0246e4d72c2eb47da3fadbd402ee))
