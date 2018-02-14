\name{vcf_example_bc}
\alias{vcf_example_bc}
\docType{data}
\title{
Data generated from VCF file with biallelic markers from a f2 backcross population
}
\description{
Simulated biallelic data set for an backcross population
}
\usage{data("vcf_example_bc")}
\format{
  An object of class \code{onemap}.
}
\details{
  A total of 142 backcross individuals were genotyped with 25 markers. The data
  was generated from a VCF file. It contains chromossome and position
  informations for each marker. It is included to be used as a example in
  order to understand how to convert VCF file to OneMap input data with the functions 
  \code{vcf2raw} and \code{onemap_read_vcfR}. 
}
\examples{
data(vcf_example_bc)
plot(vcf_example_bc) 
}
\author{Cristiane Hayumi Taniguti, \email{chaytaniguti@gmail.com}}
\seealso{
  \code{\link[onemap]{read_onemap}} for details about objects of class
  \code{onemap}.
}
\keyword{datasets}
