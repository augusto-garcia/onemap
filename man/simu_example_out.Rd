\name{simu_example_out}
\docType{data}
\alias{simu_example_out}
\title{Simulated data from a outcrossing population}
\description{
  Simulated data set from a outcrossing population.
}  
\usage{data(simu_example_out)}
\format{
  The format is:
List of 11
 $ geno         : num [1:200, 1:54] 2 1 2 1 1 2 2 2 1 1 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:200] "F1_001" "F1_002" "F1_003" "F1_004" ...
  .. ..$ : chr [1:54] "M001" "M002" "M003" "M004" ...
 $ n.ind        : int 200
 $ n.mar        : int 54
 $ segr.type    : chr [1:54] "D2.16" "D2.17" "D2.17" "D1.9" ...
 $ segr.type.num: num [1:54] 7 7 7 6 1 3 3 1 7 6 ...
 $ n.phe        : int 0
 $ pheno        : NULL
 $ CHROM        : NULL
 $ POS          : NULL
 $ input        : chr "simu_example_out.raw"
 $ error        : num [1:10800, 1:4] 1.00e-05 1.00e-05 1.00e-05 1.00 3.33e-06 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:10800] "M001_F1_001" "M002_F1_001" "M003_F1_001" "M004_F1_001" ...
  .. ..$ : NULL
 - attr(*, "class")= chr [1:2] "onemap" "outcross"
}

\details{
  A simulation of a outcrossing population of 200 individuals genotyped with 54 markers. There are no missing data. There are two groups, one (Chr01) with a total of 100 cM and the other (Chr10) with 150 cM. The markers are positioned equidistant from each other.  
}
\examples{
data(simu_example_out)

# perform two-point analyses
twopts <- rf_2pts(simu_example_out)
twopts
}
\author{Cristiane Taniguti, \email{chtaniguti@usp.br}}
\seealso{
  \code{\link[onemap]{read_onemap}} and \code{\link[onemap]{read_mapmaker}}.
}
\keyword{datasets}
