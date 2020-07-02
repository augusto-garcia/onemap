#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: vcf2raw.R                                                     #
# Contains: vcf2raw                                                   #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2015, Gabriel R A Margarido                           #
#                                                                     #
# First version: 10/14/2015                                           #
# Last update: 04/04/2016                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################


##' Convert variants from a VCF file to OneMap file format
##'
##' Converts data from a standard VCF (Variant Call Format) file to the input
##' format required by OneMap, while trying to identify the appropriate marker
##' segregation patterns.
##'
##' The input VCF file must be sorted, compressed and tabix indexed. Please check
##' functions \code{bgzip} and \code{indexTabix} of package \code{Rsamtools} for
##' details.
##'
##' Each variant in the VCF file is processed independently. Only biallelic SNPs
##' and indels for diploid variant sites are considered.
##'
##' Genotype information on the parents is required for all cross types. For
##' full-sib progenies, both outbred parents must be genotyped. For backcrosses,
##' F2 intercrosses and recombinant inbred lines, the \emph{original inbred
##' lines} must be genotyped. Particularly for backcross progenies, the
##' \emph{recurrent line must be provided as the first parent} in the function
##' arguments.
##'
##' First, samples corresponding to both parents of the progeny are parsed and
##' their genotypes identified, given that their replicates are concordant above
##' a threshold given by \code{min_class}. This allows replicates of the parents
##' to be used, which is common in sequencing plates. In detail, each parent will
##' be called an heterozygote only if
##' \eqn{min\_class * number \ of \ replicates}{[min_class * number of replicates]}
##' samples or more are heterozygous. The same is valid for homozygous calls.
##' Whenever there are different genotypes among replicates, heterozygosity is
##' checked first. The default value (\code{1.0}) requires that all replicates be
##' of the same genotype. If each parent is represented by a single sample, this
##' parameter has no effect.
##'
##' Next, marker type is determined based on parental genotypes. Finally, progeny
##' genotypes are identified and output is produced. Variants for which parent
##' genotypes cannot be determined are discarded.
##'
##' Reference sequence ID and position for each variant site are stored as special
##' fields denoted \code{CHROM} and \code{POS}.
##'
##' @param input path to the input VCF file.
##' @param output path to the output OneMap file.
##' @param cross type of cross. Must be one of: \code{"outcross"} for full-sibs;
##' \code{"f2 intercross"} for an F2 intercross progeny; \code{"f2 backcross"};
##' \code{"ri self"} for recombinant inbred lines by self-mating; or
##' \code{"ri sib"} for recombinant inbred lines by sib-mating.
##' @param parent1 \code{string} or \code{vector} of \code{strings} specifying
##' sample ID(s) of the first parent.
##' @param parent2 \code{string} or \code{vector} of \code{strings} specifying
##' sample ID(s) of the second parent.
##' @param min_class a real number between 0.0 and 1.0. For each parent and each
##' variant site, defines the proportion of parent samples that must be of the
##' same genotype for it to be assigned to the corresponding parent.
##' @author Gabriel R A Margarido, \email{gramarga@@gmail.com}
##' @seealso \code{read_onemap} for a description of the OneMap file format.
##' @keywords IO
##' @examples
##'
##'   \dontrun{
##'     vcf2raw(input="your_VCF_file.vcf.gz",
##'             output="your_OneMap_file.raw",
##'             cross="your_cross_type",
##'             parent1=c("PAR1_sample1", "PAR1_sample2"),
##'             parent2=c("PAR2_sample1", "PAR2_sample2", "PAR2_sample3"),
##'             min_class=0.5) # for parent1, a single heterozygote replicate results
##'                            # in a heterozygote genotype call; for parent2, at
##'                            # least two samples have to be concordant
##'   }
##'   
##'@export

vcf2raw <- function(input = NULL, output = NULL,
                    cross = c("outcross", "f2 intercross", "f2 backcross", "ri self", "ri sib"),
                    parent1 = NULL, parent2 = NULL, min_class = 1.0) {
  .Defunct(msg = "Defunct since version 2.1.1006. See onemap_read_vcfR function.")
  if (is.null(input)) {
    stop("You must specify the input file path.")
  }
  if (!file.exists(input)) {
    stop("Input file not found.")
  }
  if (is.null(output)) {
    stop("You must specify the output file path.")
  }
  
  ## Get absolute file paths to pass on to C
  input <- normalizePath(input)
  output <- normalizePath(output, mustWork = FALSE)
  
  cross <- match.arg(cross)
  
  if (is.null(parent1) || is.null(parent2)) {
    stop("You must specify at least one sample each as parents 1 and 2.")
  }

  # convert <- .C("vcf2raw",
  #               as.character(input),
  #               as.character(output),
  #               as.character(cross),
  #               as.integer(length(parent1)),
  #               as.character(parent1),
  #               as.integer(length(parent2)),
  #               as.character(parent2),
  #               as.numeric(min_class),
  #               PACKAGE = "onemap")
}
