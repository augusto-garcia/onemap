#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# zzz.R                                                               #
# Contains: .onemapEnv                                                #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido and Marcelo Mollinari  #
# copyright (c) 2007, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 03/12/2012                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

.onemapEnv <- new.env()
assign(".map.fun",  "kosambi", envir = .onemapEnv)
# end of file

.onAttach <- function(libname, pkgname){
  msg <- paste(
    "\n\nSince version 3.2.0, updates to **OneMap** have focused primarily on maintaining accessibility and functionality.\n\n",
    
    "For more efficient use of computational resources and time, we encourage users, especially those working with large or biallelic datasets, to consider polyploid-focused packages such as MAPpoly",
    "Although designed for polyploids, MAPpoly often runs faster than OneMap when applied to diploid datasets. ",
    "In particular, we recommend using **MAPpoly instead of OneMap** in the following scenarios involving diploid species:\n\n",
    
    "- Analyses based solely on biallelic markers (e.g., SNPs).\n",
    "- Outcrossing populations such as full-sib (F1), F2, or backcross designs.\n",
    "- Datasets containing a large number of markers (>5,000).\n",
    "- Multi-population datasets (e.g., progeny derived from multiple parents; see **MAPpoly2**, and **polyOrigin**).\n\n",
    
    "**OneMap** remains the best option when working with:\n\n",
    "- Recombinant inbred line (RIL) populations.\n",
    "- Datasets containing multiallelic or dominant markers.\n\n",
    
    "If you are using OneMap (version 2.0 or later), please refer to the following publication:\n\n",
    "C. H. Taniguti, L. M. Taniguti, R. R. Amadeu, J. Lau, G. de S. Gesteira, T. de P. Oliveira, G. C. Ferreira, ",
    "G. da S. Pereira, D. Byrne, M. Mollinari, O. Riera-Lizarazu, A. A. F. Garcia, ",
    "Developing best practices for genotyping-by-sequencing analysis in the construction of linkage maps, ",
    "GigaScience, Volume 12, 2023, giad092. https://doi.org/10.1093/gigascience/giad092",
    sep = ""
  )
  
  packageStartupMessage(msg)
}

