% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{remove_inds}
\alias{remove_inds}
\title{Remove individuals from the onemap object}
\usage{
remove_inds(onemap.obj, rm.ind)
}
\arguments{
\item{onemap.obj}{object of class onemap}

\item{rm.ind}{vector of characters with individuals names}
}
\value{
An object of class \code{onemap} without the selected individuals, 
i.e., a list with the following
components: \item{geno}{a matrix with integers indicating the genotypes
read for each marker. Each column contains data for a marker and each row
represents an individual.} \item{n.ind}{number of individuals.}
\item{n.mar}{number of markers.} \item{segr.type}{a vector with the
segregation type of each marker, as \code{strings}.} \item{segr.type.num}{a
vector with the segregation type of each marker, represented in a
simplified manner as integers, i.e. 1 corresponds to markers of type
\code{"A"}; 2 corresponds to markers of type \code{"B1.5"}; 3 corresponds
to markers of type \code{"B2.6"}; 4 corresponds to markers of type
\code{"B3.7"}; 5 corresponds to markers of type \code{"C.8"}; 6 corresponds
to markers of type \code{"D1"} and 7 corresponds to markers of type
\code{"D2"}. Markers for F2 intercrosses are coded as 1; all other crosses
are left as \code{NA}.} \item{input}{the name of the input file.}
\item{n.phe}{number of phenotypes.} \item{pheno}{a matrix with phenotypic
values. Each column contains data for a trait and each row represents an
individual.}
}
\description{
Remove individuals from the onemap object
}
\author{
Cristiane Taniguti, \email{chtaniguti@tamu.edu}
}
