\name{rf.graph.table}
\alias{rf.graph.table}
\title{Plots pairwise recombination fractions and LOD Scores using a color scale.}
\description{
  Plots a matrix of pairwise recombination fractions (under the
  diagonal) and LOD Scores (upper the diagonal) using a color scale. Any
  value of the matrix can be easily accessed using an interactive Tcl-Tk
  interface, helping the user to check for possible problems.
}
\usage{
rf.graph.table(input.seq, scale=1, axis.cex=1, main, inter=TRUE)
}
\arguments{
  \item{input.seq}{an object of class \code{sequence} with a predefined order.}
  \item{scale}{controls the plot size. If \code{inter == FALSE} this value is not used.}
  \item{axis.cex}{the magnification to be used for axis annotation.}
  \item{main}{the title for no interactive plot, i.e. it is only
    used if \code{inter == FALSE}.}
  \item{inter}{logical. If \code{TRUE}, an interactive graphic is plotted. Otherwise, a default graphic device is used.}
}
\details{
  The color scale varies from red (small distances or big LODs) to
  dark blue. When clicking on a cell, a dialog box is displayed with some
  information about corresponding markers for that cell (line
  \eqn{\times} column). The informations are: \eqn{i}) the name of the markers;
  \eqn{ii}) the number of the markers on the data set; \eqn{iii}) the segregation
  types; \eqn{iv}) the recombination fraction between the  markers and
  \eqn{v}) the LOD-Score for each possible linkage phase calculated via
  two-point analysis. For neighbor
  markers, the multipoint recombination fraction is printed; otherwise,
  the two-point recombination fraction is printed. For markers of type \code{D1}
  and \code{D2}, it's impossible to calculate recombination
  fraction via two-point analysis and, therefore the corresponding cell will be
  empty. For cells on the diagonal of the matrix, the name, the number and the type of the marker
  are printed, as well as the percentage of missing data for that marker.
}

\author{Marcelo Mollinari, \email{mmollina@usp.br}}
\examples{
##outcross example
  data(example.out)
  twopt <- rf.2pts(example.out)
  all.mark <- make.seq(twopt,"all")
  groups <- group(all.mark)
  LG1 <- make.seq(groups,1)
  LG1.rcd <- rcd(LG1)
  rf.graph.table(LG1.rcd, inter=FALSE)
\dontrun{
  ##Now, using interactive Tcl-Tk
  rf.graph.table(LG1.rcd, scale=1.5, inter=TRUE)

  ##F2 example
  data(fake.f2.onemap)
  twopt <- rf.2pts(fake.f2.onemap)
  all.mark <- make.seq(twopt,"all")
  groups <- group(all.mark)

  ##"pre-allocate" an empty list of length groups$n.groups (3, in this case)
  maps.list<-vector("list", groups$n.groups)

  for(i in 1:groups$n.groups){
    ##create linkage group i  
    LG.cur <- make.seq(groups,i)
    ##ordering 
    map.cur<-order.seq(LG.cur, subset.search = "sample")
    ##assign the map of the i-th group to the maps.list 
    maps.list[[i]]<-make.seq(map.cur, "force")      
  }
  ##Plot LOD/recombination fraction matrices for each group 
  op <- par(mfrow = c(1, 3))
  for(i in 1:groups$n.groups)
    rf.graph.table(maps.list[[i]], axis.cex=.7, main=paste("Group", i),inter=FALSE)
  par(op)
}
 

}
\keyword{utilities}
