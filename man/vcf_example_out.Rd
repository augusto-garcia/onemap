\name{vcf_example_out}
\docType{data}
\alias{vcf_example_out}
\title{Data generated from VCF file with biallelic markers from a full-sib family derived from two outbred parents}
\description{
  Simulated biallelic data set for an outcross, i.e., an F1 population obtained by
  crossing two non-homozygous parents.
}
\usage{data(vcf_example_out)}
\format{
  An object of class \code{onemap}.
}
\details{
  A total of 92 F1 individuals were genotyped with 27 markers. The data
  was generated from a VCF file. It contains chromossome and position
  informations for each marker. It is included to be used as a reference in
  order to understand how to convert VCF file to OneMap input data. Also,
  it is used for the analysis in the tutorial that comes with OneMap.
}
\examples{
data(vcf_example_out)

# plot markers informations
plot(vcf_example_out)
}

\author{Cristiane Hayumi Taniguti, \email{chaytaniguti@gmail.com}}
\seealso{
  \code{\link[onemap]{read_onemap}} for details about objects of class
  \code{onemap}.
}
\keyword{datasets}
