<!-- badges: start -->
[![Build Status](https://travis-ci.com/Cristianetaniguti/onemap.svg?branch=master)](https://travis-ci.com/Cristianetaniguti/onemap) 
[![R-CMD-check](https://github.com/Cristianetaniguti/onemap/workflows/R-CMD-check/badge.svg)](https://github.com/Cristianetaniguti/onemap/actions)
[![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![codecov](https://codecov.io/gh/Cristianetaniguti/onemap/branch/master/graph/badge.svg)](https://codecov.io/gh/Cristianetaniguti/onemap)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/onemap)](https://cran.r-project.org/package=onemap)
[![CRAN_monthly_downloads](https://cranlogs.r-pkg.org/badges/onemap)](https://cranlogs.r-pkg.org/badges/onemap)
  <!-- badges: end -->
  
# OneMap <img src="https://user-images.githubusercontent.com/7572527/119237022-0b19a400-bb11-11eb-9d45-228a59f22a1a.png" align="right" width="200"/>

`OneMap` is a software for constructing genetic maps in experimental
crosses: full-sib, RILs, F2, and backcrosses. It was developed by
Gabriel R A Margarido, Marcelo Mollinari and A Augusto F Garcia. Later on, 
Rodrigo R Amadeu, Cristiane H Taniguti, and Getulio C. Ferreira joined the project.

It has been available on CRAN for several years
(https://cran.r-project.org/package=onemap). Its last version was
updated on 2020-02-17. CRAN has OneMap's stable version, which is
recommended for most users.

This GitHub page has its version under development. New functions will
be added (experimental work) and, once it is done, we will synchronize
the repositories and add them to CRAN.

We worked very hard to release a new stable version allowing users to
analyze data sets with markers based on sequencing technologies, such
as Illumina, GBS, etc.

`OneMap` comprises a set of functions that allow users to build a
linkage map. Some functions are used internally by the package, and
should not be used directly.

# How to install

## From CRAN (stable version)

It is easy, just type (within R):

```R
setRepositories(ind = 1:2)
install.packages("onemap", dependencies=TRUE)
```

You also can use the console menus: _Packages -> Install
package(s)_. After clicking, a box will pop-up asking you to choose
the CRAN mirror. Choose the location nearest to you. Then, another box
will pop-up asking you to choose the package you want to install.
Select _onemap_ then click _OK_. The package will be
automatically installed on your computer.

`OneMap` can also be installed by downloading the appropriate files
directly at the CRAN website and following the instructions given in section `6.3 Installing Packages` of the
[R Installation and Administration](https://cran.r-project.org/doc/manuals/R-admin.pdf)
manual.

## From github (version under development)

Within R, you need to install and load the package `devtools`:

```R
install.packages("devtools")
library(devtools)
```

This will allow you to automatically build and install packages from
GitHub. If you use Windows, first install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). If you are facing problems with Rtools installation, try to do it by selecting the *Run as Administrator* option with the right mouse button. On a Mac,
you will need Xcode (available on the App Store).

Then, to install `OneMap` from GitHub (this very repo):

```R
install_github("augusto-garcia/onemap")
```

## From docker hub

`OneMap` requires several dependencies that you may not have in your system. To overcome the need of installing all of them, you can use the `OneMap` image in the docker hub. Install docker (see more about [here](https://docs.docker.com/get-started/)) and use:

```bash
docker pull cristaniguti/onemap_git:latest
```

The `OneMap` image already has the RStudio from rocker image, you can run it in your favorite browser running the following command:

```bash
docker run -p 8787:8787 -v $(pwd):/home/rstudio/ -e DISABLE_AUTH=true cristaniguti/onemap_git
```

This will make the container available in port 8787 (choose other if you prefer). The `-v` argument includes directories of your computer, in this case, the current directory (pwd) to the container. You can use `-v` several times to include several directories. After, you just need to go to your favorite browser and search for <your_localhost>:8787 (example 127.0.0.1:8787). That is it! Everything you need is there.

# Tutorials

You can read _OneMap_ tutorials going to the vignettes of the
installed package, or clicking below. Please, start with the overview,
that will guide you through other chapters.

1. [Overview](https://statgen-esalq.github.io/tutorials/onemap/Overview.html)

2. [Introduction to R](https://statgen-esalq.github.io/tutorials/onemap/Introduction_R.html)

3. [How to build a linkage map for inbred-bases populations (F2, RIL and BC)](https://statgen-esalq.github.io/tutorials/onemap/Inbred_Based_Populations.html)

4. [How to build a linkage map for outcrossing populations](https://statgen-esalq.github.io/tutorials/onemap/Outcrossing_Populations.html)

5. [A guide to build high-density linkage maps](https://cristianetaniguti.github.io/Tutorials/onemap/Quick_HighDens/High_density_maps.html)

# How to cite

Margarido, G. R. A., Souza, A. P., &38; Garcia, A. A. F. (2007). OneMap: software for genetic mapping in outcrossing species. Hereditas, 144(3), 78–79. https://doi.org/10.1111/j.2007.0018-0661.02000.x

* If you are using OneMap versions > 2.0, please cite also:

Taniguti, C. H., Taniguti, L. M., Amadeu, R. R., Mollinari, M., Da, G., Pereira, S., Riera-Lizarazu, O., Lau, J., Byrne, D., de Siqueira Gesteira, G., De, T., Oliveira, P., Ferreira, G. C., &; Franco Garcia, A. A.  Developing best practices for genotyping-by-sequencing analysis using linkage maps as benchmarks. BioRxiv. https://doi.org/10.1101/2022.11.24.517847

* If you used the HMM parallelization, please cite [BatchMap](https://github.com/bschiffthaler/BatchMap) paper too:

Schiffthaler, B., Bernhardsson, C., Ingvarsson, P. K., &; Street, N. R. (2017). BatchMap: A parallel implementation of the OneMap R package for fast computation of F1 linkage maps in outcrossing species. PLoS ONE, 12(12), 1–12. https://doi.org/10.1371/journal.pone.0189256