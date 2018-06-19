\name{onemap_example_riself}
\alias{onemap_example_riself}
\docType{data}
\title{
Simulated data from a RIL population produced by selfing.
}
\description{
Simulated biallelic data set for an \code{ri self} population.
}
\usage{data("onemap_example_riself")}
\format{
  The format is:
List of 10
 $ geno         : num [1:100, 1:68] 3 1 3 1 1 1 1 1 1 1 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:100] "ID1" "ID2" "ID3" "ID4" ...
  .. ..$ : chr [1:68] "M1" "M2" "M3" "M4" ...
 $ n.ind        : int 100
 $ n.mar        : int 68
 $ segr.type    : chr [1:68] "A.B" "A.B" "A.B" "A.B" ...
 $ segr.type.num: logi [1:68] NA NA NA NA NA NA ...
 $ n.phe        : int 0
 $ pheno        : NULL
 $ CHROM        : NULL
 $ POS          : NULL
 $ input        : chr "onemap_example_riself.raw"
 - attr(*, "class")= chr [1:2] "onemap" "riself"
}
\details{
  A total of 100 F1 individuals were genotyped for 68 markers. The data
  currently contains only genotype information (no phenotypes). It is
  included to be used as a reference in order to understand how a data
  file needs to be. 
}

\examples{
data(onemap_example_riself)
plot(onemap_example_riself)
}
\author{Cristiane Taniguti, \email{chtaniguti@usp.br}}
\seealso{
  \code{\link[onemap]{read_onemap}} for details about objects of class
  \code{onemap}.
}
\keyword{datasets}
