\name{fake_bc_onemap}
\docType{data}
\alias{fake_bc_onemap}
\title{Simulated data from a backcross population}
\description{
  Simulated data set from a backcross population.
}  
\usage{data(fake_bc_onemap)}
\format{
  An object of class \code{onemap}. 
}
\details{
  A total of 150 individuals were genotyped for 67 markers with 15\% of
  missing data. There is one quantitative phenotype to show how
  to use \code{onemap} output as \code{R\\qtl} input.
}
\examples{
data(fake_bc_onemap)

# perform two-point analyses
twopts <- rf_2pts(fake_bc_onemap)
twopts
}
\author{Marcelo Mollinari, \email{mmollina@usp.br}}
\seealso{
  \code{\link[onemap]{read_onemap}} and \code{\link[onemap]{read_mapmaker}}.
}
\keyword{datasets}
