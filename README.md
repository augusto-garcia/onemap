# OneMap

`OneMap` is a software for constructing genetic maps in experimental
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

`OneMap` comprises by a set of functions that allows users to build a
linkage map. Some functions are used internally by the package, and
should not be used directly.

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

You also can use the console menus: _Packages -> Install
package(s)_. After clicking, a box will pop-up asking you to choose
the CRAN mirror. Choose the location nearest you. Then, another box
will pop-up asking you to choose the package you want to install.
Select _onemap_ then click _OK_. The package will be
automatically installed on your computer.

`OneMap` can also be installed by downloading the appropriate files
directly at the CRAN web site and following the instructions given in
the section `6.3 Installing Packages` of the
[R Installation and Administration](http://cran.r-project.org/doc/manuals/R-admin.pdf)
manual.

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

# Tutorials

You can read _OneMap_ tutorials going to the vignettes of the
installed package, or clicking below. Please, start with the overview,
that will guide you through another chapters.

1.
[Overview](http://htmlpreview.github.com/?https://github.com/augusto-garcia/onemap/blob/master/vignettes_html/Overview.html)

2.
[Introduction to R](http://htmlpreview.github.com/?https://github.com/augusto-garcia/onemap/blob/master/vignettes_html/Introduction_R.html)

3. [Build a linkage map for inbred-bases populations (F2, RIL and BC)](http://htmlpreview.github.com/?https://raw.githubusercontent.com/augusto-garcia/onemap/master/vignettes_html/Inbred_Based_Populations.html)
