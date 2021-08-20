\name{simu_example_f2}
\docType{data}
\alias{simu_example_f2}
\title{Simulated data from a f2 intercross population}
\description{
  Simulated data set from a f2 intercross population.
}  
\usage{data(simu_example_f2)}
\format{
  The format is:
List of 11
 $ geno         : num [1:200, 1:54] 1 2 1 1 2 2 1 1 1 2 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:200] "F2_001" "F2_002" "F2_003" "F2_004" ...
  .. ..$ : chr [1:54] "M001" "M002" "M003" "M004" ...
 $ n.ind        : int 200
 $ n.mar        : int 54
 $ segr.type    : chr [1:54] "C.A" "C.A" "C.A" "C.A" ...
 $ segr.type.num: num [1:54] 7 7 7 7 4 4 7 4 4 4 ...
 $ n.phe        : int 0
 $ pheno        : NULL
 $ CHROM        : NULL
 $ POS          : NULL
 $ input        : chr "simu_example_f2.raw"
 $ error        : num [1:10800, 1:4] 1 1 1 1 1 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:10800] "M001_F2_001" "M002_F2_001" "M003_F2_001" "M004_F2_001" ...
  .. ..$ : NULL
 - attr(*, "class")= chr [1:2] "onemap" "f2"
}

\details{
  A simulation of a f2 intercross population of 200 individuals genotyped with 54 markers. There are no missing data. There are two groups, one (Chr01) with a total of 100 cM and the other (Chr10) with 150 cM. The markers are positioned equidistant from each other.  
}
\examples{
data(simu_example_f2)

# perform two-point analyses
twopts <- rf_2pts(simu_example_f2)
twopts
}
\author{Cristiane Taniguti, \email{chtaniguti@usp.br}}
\seealso{
  \code{\link[onemap]{read_onemap}} and \code{\link[onemap]{read_mapmaker}}.
}
\keyword{datasets}
