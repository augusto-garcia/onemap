\name{onemap_example_f2}
\alias{onemap_example_f2}
\docType{data}
\title{
Simulated data from a F2 population
}
\description{
 Simulated data set from a F2 population.
}
\usage{data("onemap_example_f2")}
\format{
  The format is:
List of 10
 $ geno         : num [1:200, 1:66] 1 3 2 2 1 0 3 1 1 3 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:200] "IND1" "IND2" "IND3" "IND4" ...
  .. ..$ : chr [1:66] "M1" "M2" "M3" "M4" ...
 $ n.ind        : int 200
 $ n.mar        : int 66
 $ segr.type    : chr [1:66] "A.H.B" "C.A" "D.B" "C.A" ...
 $ segr.type.num: num [1:66] 1 3 2 3 3 2 1 3 2 1 ...
 $ n.phe        : int 1
 $ pheno        : num [1:200, 1] 37.6 36.4 37.2 35.8 37.1 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr "Trait_1"
 $ CHROM        : NULL
 $ POS          : NULL
 $ input        : chr "/home/cristiane/R/x86_64-pc-linux-gnu-library/3.4/onemap/extdata/onemap_example_f2.raw"
 - attr(*, "class")= chr [1:2] "onemap" "f2"
}
\details{
  A total of 200 individuals were genotyped for 66 markers (36
  co-dominant, i.e. a, ab or b and 30 dominant i.e. c or a and d or b) with 15\% of missing data. 
  There is one quantitative phenotype to show how
  to use \code{onemap} output as \code{R\\qtl} and \code{QTL Cartographer} input. Also, it is used
  for the analysis in the tutorial that comes with OneMap.
}

\examples{
data(onemap_example_f2)
plot(onemap_example_f2)
}
\keyword{datasets}
