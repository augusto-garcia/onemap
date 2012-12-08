\name{onemap-package}
\alias{onemap-package}
\alias{onemap}
\docType{package}
\title{
  Software for constructing genetic maps
}
\description{
  In its earlier version, the software implemented the methodology
  proposed by \cite{Wu et al.} (2002a), which
  uses the maximum likelihood approach to estimate linkage and linkage
  phases for a mixed set of different marker types, to building genetic
  maps in outcrossing species.
  After, it was modified to also using a Hidden Markov Model  approach
  for constructing multipoint maximum likelihood linkage maps (Wu et
  al. 2002b). It was then applied in several studies, such as Garcia et
  al. (2006), Oliveira et al. (2007) and Oliveira et
  al. (2008). Nowadays, the latest versions (2.0-3) was  fully modified
  to also handle with backcrosses, F_2 and recombinant inbred lines
  (RIL) populations, allowing fully integration with software for QTL
  mapping, such as R/qtl and QTL Cartographer. OneMap can read and
  export files to this packages, and also from the widely used software
  MAPMAKER.
}

\details{
  \tabular{ll}{
    Package: \tab onemap\cr
    Type: \tab Package\cr
    Version: \tab 2.0-3\cr
    Date: \tab 2012-12-06\cr
    License: \tab GNU GENERAL PUBLIC LICENSE (Version 3)\cr
  }
  
  Usage of \pkg{onemap} is completely described in a tutorial
  distributed with the package. You can look for it in directory
  \code{/doc} of the package distribution.

  The most important functions are:
  \itemize{
    \item \code{\link[onemap]{read.outcross}} for data importing.

    \item \code{\link[onemap]{read.mapmaker}} for importing data from a
    MAPMAKER raw file.

    \item \code{\link[onemap]{rf.2pts}} to perform two-point
    analyses between all pairs of markers.

    \item \code{\link[onemap]{marker.type}} to check the segregation type of
    a marker.

    \item \code{\link[onemap]{make.seq}} to define a sequence of markers,
    which are the input of most mapping functions.
     
    \item \code{\link[onemap]{group}} to assign markers to linkage
    groups.

    \item \code{\link[onemap]{compare}} to compare all possible orders
    of markers in a sequence.

    \item \code{\link[onemap]{try.seq}} to try a given marker in every
    position of a framework of mapped markers.

    \item \code{\link[onemap]{order.seq}} to automate the process of mapping
    a sequence of markers, combining compare and try.seq functions.

    \item \code{\link[onemap]{ripple.seq}} to check a set of mapped markers,
    looking for (plausible) alternative orders.

    \item \code{\link[onemap]{map}} to construct a map for a sequence
    in a given order.
  }
}
\note{
  See the complete tutorial distributed along with this package for
  complete examples (located in directory /doc of the package
  distribution).
}
\author{Gabriel R A Margarido, \email{gramarga@gmail.com},
  Marcelo Mollinari, \email{mmollina@usp.br} and
  Antonio Augusto F Garcia, \email{augusto.garcia@usp.br}}
\references{
  Basten, C. J., Weir, B. S. and Zeng, Z.-B. (2005) \emph{QTL
  Cartographer Version 1.17: A Reference Manual and Tutorial for QTL
  Mapping}.
  
  Broman, K. W., Wu, H., Churchill, G., Sen, S., Yandell, B. (2008)
  \emph{qtl: Tools for analyzing QTL experiments} R package version
  1.09-43

  Buetow KH, Chakravarti A (1987). Multipoint gene mapping using
  seriation. I. General methods. \emph{Am J Hum} Genet 41: 180-188. 
  
  Doerge, R. W. (1996) Constructing genetic maps by rapid chain
  delineation. \emph{Journal of Quantitative Trait Loci} 2: 121-132.
  
  Garcia, A. A. F., Kido, E. A., Meza, A. N., Souza, H. M. B., Pinto,
  L. R., Pastina, M. M., Leite, C. S., Silva, J. A. G., Ulian, E. C.,
  Figueira, A. V. O. and Souza, A. P. (2006) Development of an integrated
  genetic map of a sugarcane (\emph{Saccharum} spp.) commercial cross,
  based on a maximum-likelihood approach for estimation of linkage and
  linkage phases. \emph{Theoretical and Applied Genetics} 112: 298-314.
  
  Haldane, J. B. S. (1919) The combination of linkage values and the
  calculation of distance between the loci of linked
  factors. \emph{Journal of Genetics} 8: 299-309.
  
  Jiang, C. and Zeng, Z.-B. (1997). Mapping quantitative trait loci with
  dominant and missing markers in various crosses from two inbred lines.
  \emph{Genetica} 101: 47-58.

  Kosambi, D. D. (1944) The estimation of map distance from
  recombination values. \emph{Annuaire of Eugenetics} 12: 172-175.
  
  Lander, E. S. and Green, P. (1987). Construction of multilocus genetic
  linkage maps in humans. \emph{Proc. Natl. Acad. Sci. USA} 84: 2363-2367.

  Lander, E. S., Green, P., Abrahamson, J., Barlow, A., Daly, M. J., Lincoln,
  S. E. and Newburg, L. (1987) MAPMAKER: An interactive computer package for
  constructing primary genetic linkage maps of experimental and natural
  populations. \emph{Genomics} 1: 174-181.

  Lincoln, S. E., Daly, M. J. and Lander, E. S. (1993) Constructing
  genetic linkage maps with MAPMAKER/EXP Version 3.0: a tutorial and
  reference manual. \emph{A Whitehead Institute for Biomedical Research
    Technical Report}.
  
  Margarido, G. R. A., Souza, A.P. and Garcia, A. A. F. (2007) OneMap:
  software for genetic mapping in outcrossing species. \emph{Hereditas}
  144: 78-79.

  Mollinari, M., Margarido, G. R. A., Vencovsky, R. and Garcia,
  A. A. F. (2009) Evaluation of algorithms used to order markers on
  genetic maps. \emph{Heredity} 103: 494-502
  
  Oliveira, K. M., Pinto, L. R., Marconi, T. G., Margarido, G. R. A.,
  Pastina, M. M., Teixeira, L. H. M., Figueira, A. V. O., Ulian, E. C.,
  Garcia, A. A. F. and Souza, A. P. (2007) Functional genetic linkage map
  based on EST-markers for a sugarcane (\emph{Saccharum} spp.)
  commercial cross. \emph{Molecular Breeding} 20: 189-208.

  Oliveira, E. J., Vieira, M. L. C., Garcia, A. A. F., Munhoz,
  C. F.,Margarido, G. R.A., Consoli, L., Matta, F. P., Moraes, M. C.,
  Zucchi, M. I., and Fungaro,M. H. P. (2008) An Integrated Molecular Map
  of Yellow Passion Fruit Based on Simultaneous Maximum-likelihood
  Estimation of Linkage and Linkage Phases
  \emph{J. Amer. Soc. Hort. Sci.} 133: 35-41.

  Tan, Y., Fu, Y. (2006). A novel method for estimating linkage
  maps. \emph{Genetics} 173: 2383-2390.

  Van Os H, Stam P, Visser RGF, Van Eck HJ (2005). RECORD: a novel
  method for ordering loci on a genetic linkage map. \emph{Theor Appl Genet}
  112: 30-40.
  
  Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002a) Simultaneous
  maximum likelihood estimation of linkage and linkage phases in
  outcrossing species. \emph{Theoretical Population Biology} 61:
  349-363.

  Wu, R., Ma, C.-X., Wu, S. S. and Zeng, Z.-B. (2002b). Linkage
  mapping of sex-specific differences. \emph{Genetical Research}
  79: 85-96
}
\keyword{package}

