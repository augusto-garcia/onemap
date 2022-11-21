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

# onemap 2.2.0

* split_onemap function to split onemap objects
* filter_missing function to filter markers in onemap object by missing data
* map parallelized with Batchmap method
* Option rm_unlinked in map to remove markers which not reach the linkage criterias when estimating genetic distances by HMM
* All ordering algorithms now presents the possibility to parallelize the last step (map) and remove automatically unlinked markers and repeat the ordering procedure
* suggest_lod function now works also with sequence objects
* Remove argument hmm from mds_onemap
* Solve numerical problems with small negative LOD values in rf_2pts 
* Update all RDatas to new format and compressed form
* Function map_avoid_unlinked to automatically remove the unlinked markers and restart the map
* create_depths_profile plots allele counts and genotypes from vcf, onemap and errors
* runpedsim  wrapper function for PedigreeSim software
* pedsim2raw converts PedigreeSim outputs to onemap raw file
* pedsim2vcf converts PedigreeSim outputs to VCF file simulating counts with updog or negative binomial
* Replace all get() by nothing
* Export parents and progeny haplotypes with write_haplotypes.R
* Draw graphics for progeny haplotypes
* F2 intercross populations now uses the same HMM than outcrossing populations. Modification needed to correct infer dominant markers and also to export progeny haplotypes.
* Function vcf2progeny_haplotypes to convert phased VCF in onemap_progeny_haplotypes object. It make possible to draw the haplotypes for phased VCFs.
* Function seq_by_type split sequence by marker type
* Internal function map_save_ram creates new onemap object for only the evaluated sequence and, after finished process that require too much RAM memory, it return the complete onemap object to user. This save RAM memory in parallelization process.
* Vignettes update
* Documentation update

# onemap 2.3.0

* Remove functions to perform simulations
* Remove updog and PedigreeSim dependency
* Edit haplotypes graphics

# onemap 2.4.0

* Add testthat tests
* Bug fix in functions

# onemap 2.5.0

* Add group_upgma function

# onemap 2.6.0

* HMM parallelization also available for Windows systems
* Add simulated data for tests purpose

# onemap 2.6.5

* Removing MDSMap from dependencies

# onemap 2.7.0

* Updated vignettes
* Github workflow
* Include hmm=FALSE option for ordering algorithms

# onemap 2.8.0

* Updated vignettes
* Add parallelization.type argument to user choose between PSOCK and FORK
* vcfR as import
* remove unlist of rf_2pts

# onemap 2.8.2

* Replacing is by inherits to fix issues in R dev versions
