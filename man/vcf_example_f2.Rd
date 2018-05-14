\name{vcf_example_f2}
\alias{vcf_example_f2}
\title{Data generated from VCF file with biallelic markers from a f2 intercross population}
\description{
  Simulated biallelic data set for an f2 population
}
\usage{data(vcf_example_f2)}
\format{
  An object of class \code{onemap}.
}
\details{
  A total of 192 F2 individuals were genotyped with 25 markers. The data
  was generated from a VCF file. It contains chromossome and position
  informations for each marker. It is included to be used as a reference in
  order to understand how to convert VCF file to OneMap input data. Also,
  it is used for the analysis in the tutorial that comes with OneMap.
}
\examples{
data(vcf_example_f2)

# plot markers informations
plot(vcf_example_f2)
}

\author{Cristiane Hayumi Taniguti, \email{chaytaniguti@gmail.com}}
\seealso{
  \code{\link[onemap]{read_onemap}} for details about objects of class
  \code{onemap}.
}
\keyword{datasets}
