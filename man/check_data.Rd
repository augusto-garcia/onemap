% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_data}
\alias{check_data}
\title{Onemap object sanity check}
\usage{
check_data(x)
}
\arguments{
\item{x}{an object of class \code{onemap}}
}
\value{
if consistent, returns 0. If not consistent, returns a 
        vector with a number of tests, where \code{TRUE} indicates
        a failed test.
}
\description{
Based on MAPpoly check_data_sanity function by Marcelo Mollinari
}
\examples{

data(onemap_example_bc)
check_data(onemap_example_bc)


}
\author{
Cristiane Taniguti, \email{chtaniguti@tamu.edu}
}
