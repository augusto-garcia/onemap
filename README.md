# OneMap

_OneMap_ is a software for constructing genetic maps in experimental
crosses: full-sib, RILs, F2 and backcrosses. It was developed by
Gabriel R A Margarido, Marcelo Mollinari and A Augusto F Garcia.

It has been available on CRAN from several years
(http://cran.r-project.org/web/packages/onemap/index.html). Its last
version was updated on 2013-09-09. CRAN has OneMap's stable version,
which is recommend for most users.

This github page has his version under development, and is a fork of
Marcelo's repository. We do not intend to have a real fork, in the
sense that there will be no two different versions maintained by
different groups. Only for my convenience, I will try to add some new
functions here (experimental work) and, once it is done, we will
synchronize the repositories and add it to CRAN.

# How to install

## From CRAN (stable version)

It is easy, just type (within R):

```R
install.packages("onemap", dependencies=TRUE)
```

Some Linux users reported the error message below:

```R
ERROR: dependency ‘tkrplot’ is not available for package ‘onemap’
```

To fix it, in a terminal (outside R), install `r-cran-tkrplot`:

```R
sudo apt-get install r-cran-tkrplot
```

Then, go back to `R` and install `OneMap` as mentioned above.

## From github (version under development)

Within R, you need to install and load the package `devtools`:

```R
install.packages("devtools")
library(devtools)
```

This will allow you to automatically build and install packages from
github. To install `OneMap` from this repo:

```R
install_github("augusto-garcia/onemap")
```


