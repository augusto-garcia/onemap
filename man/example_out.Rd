\name{example_out}
\docType{data}
\alias{example_out}
\title{Data from a full-sib family derived from two outbred parents}
\description{
  Simulated data set for an outcross, i.e., an F1 population obtained by
  crossing two non-homozygous parents.
}
\usage{data(example_out)}
\format{
  An object of class \code{onemap}.
}
\details{
  A total of 100 F1 individuals were genotyped for 30 markers. The data
  currently contains only genotype information (no phenotypes). It is
  included to be used as a reference in order to understand how a data
  file needs to be. Also, it is used for the analysis in the tutorial
  that comes with OneMap.
}
\examples{
data(example_out)

# perform two-point analyses
twopts <- rf_2pts(example_out)
twopts
}
\author{Gabriel R A Margarido, \email{gramarga@gmail.com}}
\seealso{
  \code{\link[onemap]{read_onemap}} for details about objects of class
  \code{onemap}.
}
\keyword{datasets}
