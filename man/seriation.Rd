% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seriation.R
\name{seriation}
\alias{seriation}
\title{Seriation}
\usage{
seriation(
  input.seq,
  LOD = 0,
  max.rf = 0.5,
  tol = 1e-04,
  rm_unlinked = TRUE,
  size = NULL,
  overlap = NULL,
  phase_cores = 1,
  hmm = TRUE,
  parallelization.type = "PSOCK",
  verbose = TRUE
)
}
\arguments{
\item{input.seq}{an object of class \code{sequence}.}

\item{LOD}{minimum LOD-Score threshold used when constructing the pairwise
recombination fraction matrix.}

\item{max.rf}{maximum recombination fraction threshold used as the LOD
value above.}

\item{tol}{tolerance for the C routine, i.e., the value used to evaluate
convergence.}

\item{rm_unlinked}{When some pair of markers do not follow the linkage criteria, 
if \code{TRUE} one of the markers is removed and ug is performed again.}

\item{size}{The center size around which an optimum is to be searched}

\item{overlap}{The desired overlap between batches}

\item{phase_cores}{The number of parallel processes to use when estimating
the phase of a marker. (Should be no more than 4)}

\item{hmm}{logical defining if the HMM must be applied to estimate multipoint
genetic distances}

\item{parallelization.type}{one of the supported cluster types. This should 
be either PSOCK (default) or FORK.}

\item{verbose}{A logical, if TRUE it output progress status
information.}
}
\value{
An object of class \code{sequence}, which is a list containing the
following components: \item{seq.num}{a \code{vector} containing the
(ordered) indices of markers in the sequence, according to the input file.}
\item{seq.phases}{a \code{vector} with the linkage phases between markers
in the sequence, in corresponding positions. \code{-1} means that there are
no defined linkage phases.} \item{seq.rf}{a \code{vector} with the
recombination frequencies between markers in the sequence. \code{-1} means
that there are no estimated recombination frequencies.}
\item{seq.like}{log-likelihood of the corresponding linkage map.}
\item{data.name}{name of the object of class \code{onemap} with the raw
data.} \item{twopt}{name of the object of class \code{rf_2pts} with the
2-point analyses.}
}
\description{
Implements the marker ordering algorithm \emph{Seriation} (\cite{Buetow &
Chakravarti, 1987}).
}
\details{
\emph{Seriation} is an algorithm for marker ordering in linkage groups. It
is not an exhaustive search method and, therefore, is not computationally
intensive. However, it does not guarantee that the best order is always
found. The only requirement is a matrix with recombination fractions
between markers.

NOTE: When there are to many pairs of markers with the same value in the
recombination fraction matrix, it can result in ties during the ordination
process and the \emph{Seriation} algorithm may not work properly. This is
particularly relevant for outcrossing populations with mixture of markers
of type \code{D1} and \code{D2}. When this occurs, the function shows the
following error message: \code{There are too many ties in the ordination
process - please, consider using another ordering algorithm}.

After determining the order with \emph{Seriation}, the final map is
constructed using the multipoint approach (function
\code{\link[onemap]{map}}).
}
\examples{

\donttest{
  ##outcross example
  data(onemap_example_out)
  twopt <- rf_2pts(onemap_example_out)
  all_mark <- make_seq(twopt,"all")
  groups <- group(all_mark)
  LG3 <- make_seq(groups,3)
  LG3.ser <- seriation(LG3)
}
}
\references{
Buetow, K. H. and Chakravarti, A. (1987) Multipoint gene
mapping using seriation. I. General methods. \emph{American Journal of
Human Genetics} 41: 180-188.

Mollinari, M., Margarido, G. R. A., Vencovsky, R. and Garcia, A. A. F.
(2009) Evaluation of algorithms used to order markers on genetics maps.
\emph{Heredity} 103: 494-502.
}
\seealso{
\code{\link[onemap]{make_seq}}, \code{\link[onemap]{map}}
}
\author{
Gabriel R A Margarido, \email{gramarga@gmail.com}
}
\keyword{utilities}
