\name{mapmaker_example_bc}
\docType{data}
\alias{mapmaker_example_bc}
\title{Simulated data from a backcross population}
\description{
  Simulated data set from a backcross population.
}  
\usage{data(mapmaker_example_bc)}
\format{
  An object of class \code{onemap}. 
}
\details{
  A total of 150 individuals were genotyped for 67 markers with 15\% of
  missing data. There is one quantitative phenotype to show how
  to use \code{onemap} output as \code{R\\qtl} input.
}
\examples{
data(mapmaker_example_bc)

# perform two-point analyses
twopts <- rf_2pts(mapmaker_example_bc)
twopts
}
\author{Marcelo Mollinari, \email{mmollina@usp.br}}
\seealso{
  \code{\link[onemap]{read_onemap}} and \code{\link[onemap]{read_mapmaker}}.
}
\keyword{datasets}
