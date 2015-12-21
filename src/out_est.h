/*
  OneMap: software for genetic mapping in outcrossing species
  Copyright (C) 2007-2015 Gabriel R A Margarido and Marcelo Mollinari

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
  File: out_est.h 
  
  Description: Header for a set of functions to compute the
  recombination fraction in outcross experimental populations. These
  functions contain the EM algorithms for all possible combination of
  types of markers (A, B1, B2, B3, C, D1 and D2). For more detail
  refer to Wu 2002.

  Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous
  maximum likelihood estimation of linkage and linkage phases in
  outcrossing species. Theoretical Population Biology 61: 349-363.

  Written by Marcelo Mollinari

  Escola Superior de Agricultura "Luiz de Queiroz"
  Departamento de Genética - São Paulo, Brazil
  Contact: mmollina@usp.br
  First version: 09/2015
  Last update: 10/2015
*/

#include <Rcpp.h>
#include <R_ext/PrtUtil.h>
#include "utils.h"
using namespace Rcpp;
using namespace std;
#define TOL 1e-05
#define LN3 1.098612288668109
#define LN4 1.38629436111989
#define LN_75 -0.28768207245178

Rcpp::NumericVector rf_A_A(Rcpp::NumericMatrix n,
			     int n_ind,
			   int mis);
Rcpp::NumericVector rf_A_B1(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis);
Rcpp::NumericVector rf_A_B2(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis);
Rcpp::NumericVector rf_A_B3(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis);
Rcpp::NumericVector rf_A_C(Rcpp::NumericMatrix n,
			   int n_ind,
			   int mis);
Rcpp::NumericVector rf_A_D1(Rcpp::NumericMatrix n,
			   int n_ind,
			    int mis);
Rcpp::NumericVector rf_A_D2(Rcpp::NumericMatrix n,
			   int n_ind,
			    int mis);
Rcpp::NumericVector rf_B1_B1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B1_B2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B1_B3(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B1_C(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B1_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B1_D2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B2_B2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B2_B3(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B2_C(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis);
Rcpp::NumericVector rf_B2_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B2_D2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B3_B3(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B3_C(Rcpp::NumericMatrix n,
			     int n_ind,
			    int mis);
Rcpp::NumericVector rf_B3_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_B3_D2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_C_C(Rcpp::NumericMatrix n,
			     int n_ind,
			   int mis);
Rcpp::NumericVector rf_C_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			    int mis);
Rcpp::NumericVector rf_C_D2(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis);
Rcpp::NumericVector rf_D1_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
Rcpp::NumericVector rf_D2_D2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis);
