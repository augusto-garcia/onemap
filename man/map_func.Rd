\name{map_func}
\alias{map_func}
\alias{haldane}
\alias{kosambi}
\title{Mapping functions Haldane and Kosambi}
\description{
  Functions to convert recombination fractions to distance in \code{cM} (centiMorgans).
}
\usage{
haldane(rcmb)
kosambi(rcmb)
}
\arguments{
  \item{rcmb}{A recombination fraction between two markers, i.e., a number between
                \code{0} and \code{0.5}.}
}
\details{
  Haldane mapping function is defined as
  \deqn{d_{M} = -\frac{1}{2}\ln(1-2r),}{dM = -(1/2) ln(1-2r),}
  for \eqn{0 \leq r \leq 0.5}, where \code{r} stands for the recombination fraction
  in \code{rcmb}. Kosambi mapping function is
  \deqn{d_{M} = \frac{1}{4}\ln \left[ \frac{1+2r}{1-2r} \right],}{dM =
    (1/4) ln [(1+2r)/(1-2r),} for \eqn{0 \leq r \leq 0.5}, where
  \code{r} is defined as above.
}
\value{
  Both functions return a number with a distance measured in \code{cM}.
}
\references{
  Haldane, J. B. S. (1919) The combination of linkage values and the
  calculation of distance between the loci of linked
  factors. \emph{Journal of Genetics} 8: 299-309.

  Kosambi, D. D. (1944) The estimation of map distance from
  recombination values. \emph{Annuaire of Eugenetics} 12: 172-175.
}
\author{Gabriel R A Margarido, \email{gramarga@gmail.com}}
%\seealso{ }
\examples{
  # little difference for small recombination fractions
  haldane(0.05)
  kosambi(0.05)

  # greater difference as recombination fraction increases
  haldane(0.35)
  kosambi(0.35)
}
\keyword{arith}
