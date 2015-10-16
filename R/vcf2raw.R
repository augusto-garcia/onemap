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
# Last update: 10/14/2015                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################


##' Convert variants from a VCF file to OneMap file format
##' 
##' Converts data from a standard VCF (Variant Call Format) file to the input
##' format required by OneMap, while trying to identify the appropriate
##' segregation pattern.
##'
##' Each variant in the VCF file is processed independently. Only biallelic SNPs
##' and indels for diploid variant sites are considered.
##'
##' First, samples corresponding to both parents of the progeny are parsed and
##' their genotypes identified, given that they are concordant above a provided
##' threshold. Marker type (B3, D1 or D2) is determined based on parental genotypes.
##' Finally, progeny genotypes are identified and output is produced. Variants
##' for which parent genotypes cannot be determined are discarded.
##'
##' Reference sequence ID and position for each variant site are stored as special
##' phenotypes denoted \code{CHROM} and \code{POS}.
##' 
##' @param input path to the input VCF file.
##' @param output path to the output OneMap file.
##' @param parent1 \code{string} or \code{vector} of \code{strings} specifying
##' sample ID(s) of the first progeny parent.
##' @param parent2 \code{string} or \code{vector} of \code{strings} specifying
##' sample ID(s) of the second progeny parent.
##' @param min_class a real number between 0.0 and 1.0. For each parent and each
##' variant site, at least \eqn{parent sample count \times min_class} of the parent
##' samples must be of the same genotype for it to be assigned to the corresponding
##' parent.
##' @author Gabriel R A Margarido, \email{gramarga@@gmail.com}
##' @seealso \code{read.outcross} for a description of the OneMap file format.
##' @keywords IO
##' @examples
##' 
##'   \dontrun{
##'     vcf2raw(input="your_VCF_file.vcf", output="your_OneMap_file.raw",
##'             parent1=c("PAR1_sample1", "PAR1_sample2"),
##'             parent2=c("PAR2_sample1", "PAR2_sample2", "PAR2_sample3"),
##'             min_class=0.5)
##'   }
##'

vcf2raw <- function(input = NULL, output = NULL, parent1 = NULL, parent2 = NULL, min_class = 1.0) {
    if (is.null(input)) {
        stop("You must specify the input file path.")
    }
    if (!file.exists(input)) {
        stop("Input file not found.")
    }
    if (is.null(output)) {
        stop("You must specify the output file path.")
    }
    if (is.null(parent1) || is.null(parent2)) {
        stop("You must specify at least one sample each as parents 1 and 2.")
    }
    
    convert <- .C("vcf2raw",
                  as.character(input),
                  as.character(output),
                  as.integer(length(parent1)),
                  as.character(parent1),
                  as.integer(length(parent2)),
                  as.character(parent2),
                  as.numeric(min_class),
                  PACKAGE = "onemap")
}
