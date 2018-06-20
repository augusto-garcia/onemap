\name{onemap_example_bc}
\docType{data}
\alias{onemap_example_bc}
\title{Simulated data from a backcross population}
\description{
  Simulated data set from a backcross population.
}  
\usage{data(onemap_example_bc)}
\format{
  The format is:
List of 10
 $ geno         : num [1:150, 1:67] 1 2 1 1 2 1 2 1 1 2 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:150] "ID1" "ID2" "ID3" "ID4" ...
  .. ..$ : chr [1:67] "M1" "M2" "M3" "M4" ...
 $ n.ind        : int 150
 $ n.mar        : int 67
 $ segr.type    : chr [1:67] "A.H" "A.H" "A.H" "A.H" ...
 $ segr.type.num: logi [1:67] NA NA NA NA NA NA ...
 $ n.phe        : int 1
 $ pheno        : num [1:150, 1] 40.8 39.5 37.9 34.2 38.9 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr "Trait_1"
 $ CHROM        : NULL
 $ POS          : NULL
 $ input        : chr "onemap_example_bc.raw"
 - attr(*, "class")= chr [1:2] "onemap" "backcross"
}

\details{
  A total of 150 individuals were genotyped for 67 markers with 15\% of
  missing data. There is one quantitative phenotype to show how
  to use \code{onemap} output as \code{R\\qtl} input.
}
\examples{
data(onemap_example_bc)

# perform two-point analyses
twopts <- rf_2pts(onemap_example_bc)
twopts
}
\author{Marcelo Mollinari, \email{mmollina@usp.br}}
\seealso{
  \code{\link[onemap]{read_onemap}} and \code{\link[onemap]{read_mapmaker}}.
}
\keyword{datasets}
