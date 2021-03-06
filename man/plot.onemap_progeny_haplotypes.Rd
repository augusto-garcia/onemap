% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_haplotypes.R
\name{plot.onemap_progeny_haplotypes}
\alias{plot.onemap_progeny_haplotypes}
\title{Plots progeny haplotypes}
\usage{
\method{plot}{onemap_progeny_haplotypes}(
  x,
  col = NULL,
  position = "stack",
  show_markers = TRUE,
  main = "Genotypes",
  ncol = 4,
  ...
)
}
\arguments{
\item{x}{object of class onemap_progeny_haplotypes}

\item{col}{Color of parentes' homologous.}

\item{position}{"split" or "stack"; if "split" (default) the parents' homologous are plotted separately. if "stack" the parents' homologous are plotted together.}

\item{show_markers}{logical; if  \code{TRUE}, the markers (default) are plotted.}

\item{main}{An overall title for the plot; default is \code{NULL}.}

\item{ncol}{number of columns of the facet_wrap}

\item{...}{currently ignored}
}
\description{
Plots progeny haplotypes
}
\author{
Getulio Caixeta Ferreira, \email{getulio.caifer@gmail.com}

Cristiane Taniguti, \email{chtaniguti@usp.br}
}
