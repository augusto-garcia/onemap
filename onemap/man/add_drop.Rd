\name{add_drop}
\alias{add_drop}
\alias{add.marker}
\alias{drop.marker}
\title{Add or Drop Markers From a Sequence}
\description{
  Creates a new sequence by adding or dropping markers from a predetermined
  one. The markers are added in the end of the sequence.
}
\usage{
add.marker(input.seq, mrks)
drop.marker(input.seq, mrks)
}
\arguments{
\item{input.seq}{an object of class \code{sequence}.}
\item{mrks}{a vector containing the markers to be added or removed from
  the \code{sequence}.}
}

\value{
  An object of class \code{sequence}, which is a list containing the
  following components:
  \item{seq.num}{a \code{vector} containing the (ordered) indices of
    markers in the sequence, according to the input file.}
  \item{seq.phases}{a \code{vector} with the linkage phases between markers
    in the sequence, in corresponding positions. \code{-1} means that there are
    no defined linkage phases.}
  \item{seq.rf}{a \code{vector} with the recombination fractions between
    markers in the sequence. \code{-1} means that there are no estimated
    recombination fractions.}
  \item{seq.like}{log-likelihood of the corresponding linkage map.}
  \item{data.name}{name of the object of class \code{outcross} with the
    raw data.}
  \item{twopt}{name of the object of class \code{rf.2pts} with the 2-point
    analyses.}
}

\author{Marcelo Mollinari, \email{mmollina@usp.br}}
%\seealso{ }
\examples{

data(example.out)
twopt <- rf.2pts(example.out)
all.mark <- make.seq(twopt,"all")
groups <- group(all.mark)
(LG1 <- make.seq(groups,1))
(LG.aug<-add.marker(LG1, c(4,7)))
(LG.red<-drop.marker(LG1, c(1,2,3,5,6)))
 
}
\keyword{manip}
\keyword{utilities}

