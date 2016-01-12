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
  File: f2_est.h
  
  Description: Header file to a set of functions to compute the
  recombination fraction in F2 experimental populations. These
  functions contain the EM algorithms for all possible combination of
  types of markers (co-dominant and dominant).

  Written by Marcelo Mollinari

  Escola Superior de Agricultura "Luiz de Queiroz"
  Departamento de Genética - São Paulo, Brazil
  Contact: mmollina@usp.br
  First version: 09/2015
  Last update: 10/2015
*/
#include <Rcpp.h>
using namespace Rcpp;

Rcpp::NumericVector est_rf_C_C(std::vector<int> k_sub,
			       std::vector<int> k1_sub,
			       int n_ind);
Rcpp::NumericVector est_rf_C_D_43(std::vector<int> k_sub,
				  std::vector<int> k1_sub,
				  int n_ind);
Rcpp::NumericVector est_rf_C_D_51(std::vector<int> k_sub,
				  std::vector<int> k1_sub,
				  int n_ind);
Rcpp::NumericVector est_rf_D_D_43(std::vector<int> k_sub,
				  std::vector<int> k1_sub,
				  int n_ind);
Rcpp::NumericVector est_rf_D_D_51(std::vector<int> k_sub,
				  std::vector<int> k1_sub,
				  int n_ind);
Rcpp::NumericVector est_rf_D_D_43_51(std::vector<int> k_sub,
				     std::vector<int> k1_sub,
				     int n_ind);
Rcpp::NumericVector est_rf_A_A(std::vector<int> k_sub,
			       std::vector<int> k1_sub,
			       int n_ind);
