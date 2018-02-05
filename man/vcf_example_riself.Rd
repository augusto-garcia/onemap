\name{vcf_example_riself}
\alias{vcf_example_riself}
\docType{data}
\title{
Data generated from VCF file with biallelic markers from a RIL population produced by selfing
}
\description{
Simulated biallelic data set for an \code{ri self} population.
}
\usage{data("vcf_example_riself")}
\format{
  The format is:
List of 10
 $ geno         : num [1:92, 1:25] 3 3 1 3 1 3 3 1 3 1 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:92] "ID1" "ID3" "ID4" "ID5" ...
  .. ..$ : chr [1:25] "SNP16" "SNP12" "SNP17" "SNP10" ...
 $ n.ind        : int 92
 $ n.mar        : int 25
 $ segr.type    : chr [1:25] "A.B" "A.B" "A.B" "A.B" ...
 $ segr.type.num: logi [1:25] NA NA NA NA NA NA ...
 $ n.phe        : int 0
 $ pheno        : NULL
 $ CHROM        : chr [1:25] "1" "1" "1" "1" ...
 $ POS          : int [1:25] 1791 6606 9001 11326 11702 15533 17151 18637 19146 19220 ...
 $ input        : chr "vcf_example_riself.raw"
 - attr(*, "class")= chr [1:2] "onemap" "riself"
}
\details{
  A total of 92 rils individuals were genotyped with 25 markers. The data
  was generated from a VCF file. It contains chromossome and position
  informations for each marker. It is included to be used as a example in
  order to understand how to convert VCF file to OneMap input data with the functions 
  \code{vcf2raw} and \code{onemap_read_vcfR}. 
}

\examples{
data(vcf_example_riself)
plot(vcf_example_riself)
}
\author{Cristiane Hayumi Taniguti, \email{chaytaniguti@gmail.com}}
\seealso{
  \code{\link[onemap]{read_onemap}} for details about objects of class
  \code{onemap}.
}
\keyword{datasets}
