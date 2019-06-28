/*
 OneMap: software for genetic mapping in outcrossing species
 Copyright (C) 2007-9 Gabriel R A Margarido and Marcelo Mollinari

    This file is part of OneMap.

    OneMap is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    OneMap is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
*/

/*
  File: hmm_bc.cpp
  Description: Set of functions to be used with software R
               Implements the methodology of Hidden Markov Models (HMM)
	       to construct multipoint linkage maps in bc crosses

  Written by Marcelo Mollinari
  Adapted from hmm_main.c, hmm_bc.c and util.c (found in the R package qtl)
  copyright (c) 2001-10, Karl W Broman                                

  Escola Superior de Agricultura "Luiz de Queiroz"
  Departamento de Genética - São Paulo, Brazil
  Contact: mmollina@usp.br
  First version: 10/2015
  Last update: 12/2015
*/

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define THRESH 200.0


/**********************************************************************
 * 
 * est_hmm_bc
 *
 * This function re-estimates the genetic map for a chromosome in a bc cross
 *
 * geno_R       Genotype data, as a matrix. The columns represent the number
 *              of markers and the rows represent the number of individuals
 *
 * rf_R         Recombination fractions
 *  
 * verbose_R    If TRUE print tracing information.
 *
 * tol_R        Tolerance for determining convergence
 * 
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

RcppExport SEXP est_hmm_bc(SEXP geno_R, SEXP error_R, SEXP rf_R, SEXP verbose_R, SEXP tol_R);
double stepf_bc(int gen1, int gen2, double rf);
double nrecf_bc(int gen1, int gen2);
