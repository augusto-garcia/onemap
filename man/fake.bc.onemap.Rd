\name{fake.bc.onemap}
\docType{data}
\alias{fake.bc.onemap}
\title{Simulated data from a backcross population}
\description{
  Simulated data set from a backcross population.
}  
\usage{data(fake.bc.onemap)}
\format{
  An object of class \code{bc.onemap}. 
}
\details{
  A total of 150 individuals were genotyped for 67 markers with 15\% of
  missing data. There is one quantitative phenotype to show how
  to use \code{onemap} output as \code{R\\qtl} input.
}
\examples{
data(fake.bc.onemap)

# perform two-point analyses
twopts <- rf.2pts(fake.bc.onemap)
twopts
}
\author{Marcelo Mollinari, \email{mmollina@usp.br}}
\seealso{
  \code{\link[onemap]{read.mapmaker}} for details about objects of class
  \code{bc.onemap}.
}
\keyword{datasets}
