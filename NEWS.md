# onemap 2.0.5

* Organizing the files in directories to allow installation from github 
* Moving functions descriptions for roxygen style
* Solving conflit of version 2.0.4 with 2.0.5 (the one in CRAN)
* Replacing old vignettes in pdf with new ones using RMarkdown
* Adding functions to plot raw data
* Adding functions to draw graphics showing marker categories
* Adding functions to perform chi-square tests for marker segregation
* Adding functions to identify alpha under Bonferroni's argument
* Correcting documentation not compiling for RMarkdown
* Updating dependencies for other R packages

# onemap 2.0.6

* (To be included)

# onemap 2.0.7

* Updating .travis.yml for the newer standards
* Updating README in order to solve Windows user's issues.
* Updating dependencies and NAMESPACE for newer standards (R-3.4.0)
* Updating dependencies within Bioconductor packages

# onemap 2.0.8

* Fixing chi-square segregation tests for markers with missing classes

# onemap 2.1.1001

* Bug fix in function test_segregation for F2 cross
* Fix in example f2
* Codif_data changed to deal with dominant markers for F2
* Warnings removed and labels added to outcross plot
* plot.onemap palette changed
* Included C.A, D.B and A.H.B types in plot by segregation
* New argument in select_segreg
* Function create_data_bins now keeps CHROM and POS infos
* Possible to make sequences from CHROM and POS using 2pts object
* Outcrossing vcf example file are available
* New functions group_seq and print.group_seq
* Functions draw.try and draw_order are now defuncts
* Vignette for outcrossing updated
* F2 vcf example file are available
* Vignette for inbred updated

# onemap 2.1.1002

* onemap_read_vcfR function
* Examples files for rils and backcross populations 
* Vignette for inbred and outcrossing updated

# onemap 2.1.1003

* rf_graph_table updated with ggplot2 and plotly
* NAMESPACE automated with roxygen2

# onemap 2.1.1004

* write_onemap_raw function

# onemap 2.1.1005

* Vignette updated

# onemap 2.1.1006

* vcf2raw as defunct

# onemap 2.1.1007

* draw_map2 function
* Vignette update

# onemap 2.1.1008

* mds_onemap function
* Vignette update

# onemap.2.1.2

* Submitted to CRAN

# onemap 2.1.2001

* Adapt compare to rils and backcross

# onemap 2.1.3

* draw_map2 now receives data.frames
* Changing class functions to is
* Adapt try_seq to rils and backcross
* onemap_read_vcfR also accept phased VCF
* Submitted to CRAN

# onemap 3.0

* split_onemap function to split onemap objects
* onemap_read_vcfR now works for multiallelic markers
* map parallelized according with Batchmap method
* parmap function - alternative function to parallelize HMM
* HMM now uses three possible errors in emission phase: global error, genotype errors, genotype probabilities
* create_probs function makes the convertion of the three types of errors to the emission matrix
* updog_genotype function performs regenotyping with updog software
* polyRAD_genotype function performs regenotyping with polyRAD software
* create_depths_profile plots allele counts and genotypes from vcf, onemap and errors
* runpedsim function makes a interface with PedigreeSim software
* pedsim2raw converts PedigreeSim outputs to onemap raw file
* pedsim2vcf converts PedigreeSim outputs to VCF file simulating counts with updog or negative binomial
* Replace all get() by nothing
* New vignettes "Simulations" and "High Density Maps"
* version compatible with onemap_workflows and onemap_workflows_app shiny app
