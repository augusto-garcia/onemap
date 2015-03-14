\name{marker.type}
\alias{marker.type}
\title{Informs the segregation patterns of markers}
\description{
  Informs the type of segregation of all markers from an object of class
  \code{sequence}. For outcross populations it uses the notation by
  \cite{Wu et al., 2002}. For backcrosses, \eqn{F_2}{F_2} and RILs,
  it uses the traditional notation from MAPMAKER i.e. AA, AB, BB, not AA and not BB.
}
\usage{
marker.type(input.seq)
}
\arguments{
  \item{input.seq}{an object of class \code{sequence}.}
}
\details{
  The segregation types are (\cite{Wu et al., 2002}):
  \tabular{lcc}{
    Type \tab Cross \tab Segregation \cr
    A.1 \tab ab x cd \tab 1:1:1:1 \cr
    A.2 \tab ab x ac \tab 1:1:1:1 \cr
    A.3 \tab ab x co \tab 1:1:1:1 \cr
    A.4 \tab ao x bo \tab 1:1:1:1 \cr
    B1.5 \tab ab x ao \tab 1:2:1 \cr
    B2.6 \tab ao x ab \tab 1:2:1 \cr
    B3.7 \tab ab x ab \tab 1:2:1 \cr
    C8 \tab ao x ao \tab 3:1 \cr
    D1.9 \tab ab x cc \tab 1:1 \cr
    D1.10 \tab ab x aa \tab 1:1 \cr
    D1.11 \tab ab x oo \tab 1:1 \cr
    D1.12 \tab bo x aa \tab 1:1 \cr
    D1.13 \tab ao x oo \tab 1:1 \cr
    D2.14 \tab cc x ab \tab 1:1 \cr
    D2.15 \tab aa x ab \tab 1:1 \cr
    D2.16 \tab oo x ab \tab 1:1 \cr
    D2.17 \tab aa x bo \tab 1:1 \cr
    D2.18 \tab oo x ao \tab 1:1
  }
}
\value{
  Nothing is returned. Segregation types of all markers in the
  sequence are displayed on the screen.
}
\references{
  Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous
  maximum likelihood estimation of linkage and linkage phases in
  outcrossing species. \emph{Theoretical Population Biology} 61: 349-363.
}
\author{Gabriel R A Margarido, \email{gramarga@gmail.com}}
\seealso{\code{\link[onemap]{make.seq}}}
\examples{
  data(example.out)
  twopts <- rf.2pts(example.out)
  markers.ex <- make.seq(twopts,c(3,6,8,12,16,25))
  marker.type(markers.ex) # segregation type for some markers

  data(fake.f2.onemap)
  twopts <- rf.2pts(fake.f2.onemap)
  all.mrk<-make.seq(twopts, "all")
  lgs<-group(all.mrk)
  lg1<-make.seq(lgs,1)
  marker.type(lg1) # segregation type for linkage group 1

}
\keyword{manip}
\keyword{utilities}
