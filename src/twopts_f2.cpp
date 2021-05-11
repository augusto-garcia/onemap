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
  File: twopts_f2.cpp

  Description: Function to compute the recombination fraction in F2
  experimental populations.

  Written by Marcelo Mollinari

  Escola Superior de Agricultura "Luiz de Queiroz"
  Departamento de Genética - São Paulo, Brazil
  Contact: mmollina@usp.br
  First version: 09/2015
  Last update: 10/2015
*/

#include "f2_est.h"
#include <Rcpp.h>
#include "utils.h"
#include "twopts_f2.h"
using namespace Rcpp;
using namespace std;

RcppExport SEXP est_rf_f2_wrap(SEXP geno_R, SEXP mrk_R, SEXP segreg_type_R, SEXP n_ind_R, SEXP verbose_R)
{
  int n_ind = Rcpp::as<int>(n_ind_R);
  int verbose = Rcpp::as<int>(verbose_R);
  int mrk = Rcpp::as<int>(mrk_R);
  Rcpp::IntegerVector segreg_type = Rcpp::as<Rcpp::IntegerVector>(segreg_type_R);
  Rcpp::NumericVector geno = Rcpp::as<Rcpp::NumericVector>(geno_R);
  NumericMatrix z = est_rf_f2(geno, mrk, segreg_type, n_ind, verbose);
  return(wrap(z));
}

Rcpp::NumericMatrix est_rf_f2(Rcpp::NumericVector geno, 
			      int mrk, 
			      Rcpp::IntegerVector segreg_type, 
			      int n_ind, int verbose)
{
  int n_mar=((int)geno.size()/n_ind), k1, k2;
  int chunk=((n_mar*n_mar)-n_mar)/20, ct1=0, ct2=1, a1, a2, a3;
  NumericMatrix r(n_mar, n_mar);
  NumericMatrix rm(2, n_mar);
  Rcpp::NumericVector rtemp(2);
  if(mrk < 0)
    {
      if(verbose==1 && n_mar > 100)
	Rcpp::Rcout << "Computing " << ((n_mar*n_mar)-n_mar)/2 << " recombination fractions:\n\n" << "\t0%\t.";
      else if(verbose==1)
	Rcpp::Rcout << "Computing " << ((n_mar*n_mar)-n_mar)/2 << " recombination fractions ... \n";
      a1=0;
      a2=n_mar-1;
    }
  else
    {
      a1=mrk;
      a2=mrk+1;
    }
  for(int i=a1; i < a2; i++)
    {
      if(mrk < 0){
	if(verbose==1 && n_mar > 100)
	  {
	    if(ct1%(chunk/(chunk/10))==0)
	      {
		Rcpp::Rcout << ".";
		ct2++;
		if(ct2%50==0)
		  Rcpp::Rcout << "\t" << 100*ct1/(((n_mar*n_mar)-n_mar)/2) << "%\n\t" << 100*ct1/(((n_mar*n_mar)-n_mar)/2) << "%\t";
	      }
	  }
      }
      R_CheckUserInterrupt(); // check for ^C
      if(mrk < 0) a3=(i+1);
      else a3=0;
      for(int j=a3; j  < n_mar; j++)
        {
	  ct1++;
	  std::vector<int> k_sub(&geno[i*n_ind],&geno[i*n_ind+n_ind]);
	  std::vector<int> k1_sub(&geno[j*n_ind],&geno[j*n_ind+n_ind]);
	  k1=segreg_type(i); k2=segreg_type(j);
	  // Rcpp::Rcout << k1 << "--" << k2 << "\n";
	  switch(k1){
	  case 1: //C
	    switch(k2){
	    case 1:
	      rtemp = est_rf_C_C(k_sub, k1_sub, n_ind);  //Markers: Codominant - Codominant
	      break;
	    case 2:
	      rtemp = est_rf_C_D_43(k_sub, k1_sub, n_ind);  //Markers: Codominant - Dominant( (12)=4, 3 )
	      break;
	    case 3:
	      rtemp = est_rf_C_D_51(k_sub, k1_sub, n_ind);  //Markers Codominant - Dominant( (1, (23)=5)
	      break;
	    case 4:
	      rtemp = est_rf_A_A(k_sub, k1_sub, n_ind);  //Any type of markers - slower function
	      break;
	    }
	    break;
	  case 2: //D_43
	    switch(k2){
	    case 1:
	      rtemp = est_rf_C_D_43(k1_sub, k_sub, n_ind);  //Markers: Dominant( (12)=4, 3 ) - Codominant
	      break;
	    case 2:
	      rtemp = est_rf_D_D_43(k_sub, k1_sub, n_ind);  //Markers: Dominant( (12)=4, 3 ) - Dominant( (12)=4, 3 )
	      break;
	    case 3:
	      rtemp = est_rf_D_D_43_51(k_sub, k1_sub, n_ind);  //Markers: Dominant( (12)=4, 3 ) - Dominant( (1, (23)=5)
	      break;
	    case 4:
	      rtemp = est_rf_A_A(k_sub, k1_sub, n_ind);  //Any type of markers - slower function
	      break;
	    }
	    break;
	  case 3: //D51
	    switch(k2){
	    case 1:
	      rtemp = est_rf_C_D_51(k1_sub, k_sub, n_ind); //Markers  Dominant( (1, (23)=5) - Codominant
	      break;
	    case 2:
	      rtemp = est_rf_D_D_43_51(k1_sub, k_sub, n_ind); //Markers:  Dominant( (1, (23)=5) - Dominant( (12)=4, 3 )
	      break;
	    case 3:
	      rtemp = est_rf_D_D_51(k_sub, k1_sub, n_ind); //Markers:  Dominant( (1, (23)=5) - Dominant( (1, (23)=5)
	      break;
	    case 4:
	      rtemp = est_rf_A_A(k_sub, k1_sub, n_ind);  //Any type of markers - slower function
	      break;
	    }
	    break;
	  case 4:
	    switch(k2){
	    case 1: case 2: case 3: case 4:
	      rtemp = est_rf_A_A(k_sub, k1_sub, n_ind); //Any type of markers - slower function
	      break;
	    }
	    break;
	  }
	  if(mrk < 0){
	    r(j,i)=rtemp(0);
	    r(i,j)=rtemp(1);
	    rtemp(0)=rtemp(1)=0;
	  }
	  else
	    {
	      rm(0,j)=rtemp(0);
	      rm(1,j)=rtemp(1);
	      rtemp(0)=rtemp(1)=0;
	    }
	}
    }
  if(mrk < 0){
    if(verbose==1 && n_mar > 100) Rcpp::Rcout << "\t100%\n";
    return(r);
  }
  else
    {
      return(rm);
    }
}
