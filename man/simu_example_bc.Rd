\name{simu_example_bc}
\docType{data}
\alias{simu_example_bc}
\title{Simulated data from a backcross population}
\description{
  Simulated data set from a backcross population.
}  
\usage{data(simu_example_bc)}
\format{
  The format is:
List of 11
 $ geno         : num [1:200, 1:54] 1 2 1 1 2 2 2 1 1 2 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:200] "BC_001" "BC_002" "BC_003" "BC_004" ...
  .. ..$ : chr [1:54] "M001" "M002" "M003" "M004" ...
 $ n.ind        : int 200
 $ n.mar        : int 54
 $ segr.type    : chr [1:54] "A.H" "A.H" "A.H" "A.H" ...
 $ segr.type.num: num [1:54] 8 8 8 8 8 8 8 8 8 8 ...
 $ n.phe        : int 0
 $ pheno        : NULL
 $ CHROM        : NULL
 $ POS          : NULL
 $ input        : chr "sim_bc.example.raw"
 $ error        : num [1:10800, 1:2] 1 1 1 1 1 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:10800] "M001_BC_001" "M002_BC_001" "M003_BC_001" "M004_BC_001" ...
  .. ..$ : NULL
 - attr(*, "class")= chr [1:2] "onemap" "backcross"
}

\details{
  A simulation of a backcross population of 200 individuals genotyped with 54 markers. Ten percent of the genotypes are   missing data. There are two groups, one (Chr01) with a total of 100 cM and the other (Chr10) with 150 cM. The markers are positioned equidistant from each other.  
}
\examples{
data(simu_example_bc)

# perform two-point analyses
twopts <- rf_2pts(simu_example_bc)
twopts
}
\author{Cristiane Taniguti, \email{chtaniguti@usp.br}}
\seealso{
  \code{\link[onemap]{read_onemap}} and \code{\link[onemap]{read_mapmaker}}.
}
\keyword{datasets}
