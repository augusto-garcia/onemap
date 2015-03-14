\name{fake.f2.onemap}
\docType{data}
\alias{fake.f2.onemap}
\title{Simulated data from a F2 population}
\description{
  Simulated data set from a F2 population.
}  
\usage{data(fake.f2.onemap)}
\format{
  An object of class \code{f2.onemap}. 
}
\details{
  A total of 200 individuals were genotyped for 66 markers (36
  co-dominant, i.e. AA, AB or BB and 30 dominant i.e. Not AA or AA and Not
  BB or BB) with 15\% of missing data. There is one quantitative phenotype to show how
  to use \code{onemap} output as \code{R\\qtl} and \code{QTL Cartographer} input. Also, it is used
  for the analysis in the tutorial that comes with OneMap.
}
\examples{
data(fake.f2.onemap)

# perform two-point analyses
twopts <- rf.2pts(fake.f2.onemap)
twopts
}
\author{Marcelo Mollinari, \email{mmollina@usp.br}}
\seealso{
  \code{\link[onemap]{read.mapmaker}} for details about objects of class
  \code{f2.onemap}.
}
\keyword{datasets}
