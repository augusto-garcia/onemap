\name{mapmaker_example_f2}
\docType{data}
\alias{mapmaker_example_f2}
\title{Simulated data from a F2 population}
\description{
  Simulated data set from a F2 population.
}  
\usage{data(mapmaker_example_f2)}
\format{
  An object of class \code{onemap}.
}
\details{
  A total of 200 individuals were genotyped for 66 markers (36
  co-dominant, i.e. AA, AB or BB and 30 dominant i.e. Not AA or AA and Not
  BB or BB) with 15\% of missing data. There is one quantitative phenotype to show how
  to use \code{onemap} output as \code{R\\qtl} and \code{QTL Cartographer} input. Also, it is used
  for the analysis in the tutorial that comes with OneMap.
}
\examples{
data(mapmaker_example_f2)

# perform two-point analyses
twopts <- rf_2pts(mapmaker_example_f2)
twopts
}
\author{Marcelo Mollinari, \email{mmollina@usp.br}}
\seealso{
  \code{\link[onemap]{read_onemap}} and \code{\link[onemap]{read_mapmaker}}.
}
\keyword{datasets}
