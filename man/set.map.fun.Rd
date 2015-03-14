\name{set.map.fun}
\alias{set.map.fun}
\title{Defines the default mapping function}
\description{
  Defines the function that should be used to display the genetic map
  through the analysis.
}
\usage{
set.map.fun(type=c("kosambi", "haldane"))
}
\arguments{
  \item{type}{Indicates the function that should be used, which can be
    \code{"kosambi"} or \code{"haldane"}}
}

\references{
  Haldane, J. B. S. (1919) The combination of linkage values and the
  calculation of distance between the loci of linked
  factors. \emph{Journal of Genetics} 8: 299-309.

  Kosambi, D. D. (1944) The estimation of map distance from
  recombination values. \emph{Annuaire of Eugenetics} 12: 172-175.
}
\author{Marcelo Mollinari, \email{mmollina@usp.br}}
\seealso{\code{\link[onemap]{kosambi}} and \code{\link[onemap]{haldane}}}

\keyword{arith}
