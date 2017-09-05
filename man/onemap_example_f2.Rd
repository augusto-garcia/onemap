\name{onemap_example_f2}
\alias{onemap_example_f2}
\docType{data}
\title{ Simulated data from a F2 population
}
\description{
  Simulated data set from a F2 population.
}  
\usage{data(onemap_example_f2)}
\format{
  An object of class \code{onemap}.
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

# perform two-point analyses
twopts <- rf_2pts(onemap_example_f2)
twopts
}
\author{Marcelo Mollinari, \email{mmollina@usp.br}}
\seealso{
  \code{\link[onemap]{read_onemap}} and \code{\link[onemap]{read_mapmaker}}.
}
\keyword{datasets}
