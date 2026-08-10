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
  File: twopts_out.cpp 

  Description: Function to compute the recombination fraction in
  outcross experimental populations.

  Written by Marcelo Mollinari

  Escola Superior de Agricultura "Luiz de Queiroz"
  Departamento de Genética - São Paulo, Brazil
  Contact: mmollina@usp.br
  First version: 09/2015
  Last update: 11/2015
*/

#include <Rcpp.h>
#include "out_est.h"
#include "utils.h"
#include "twopts_out.h"
#define rf_TOL_min 1e-50
#define rf_TOL_max 1-1e-50
using namespace Rcpp;
using namespace std;

RcppExport SEXP est_rf_out_wrap(SEXP geno_R, SEXP mrk_R, SEXP segreg_type_R, SEXP n_ind_R, SEXP verbose_R)
{
  int n_ind = Rcpp::as<int>(n_ind_R);
  int verbose = Rcpp::as<int>(verbose_R);
  int mrk = Rcpp::as<int>(mrk_R);
  Rcpp::IntegerVector segreg_type = Rcpp::as<Rcpp::IntegerVector>(segreg_type_R);
  Rcpp::NumericVector geno = Rcpp::as<Rcpp::NumericVector>(geno_R);
  List z = est_rf_out(geno, mrk, segreg_type, n_ind, verbose);
  return(z);
}

Rcpp::List est_rf_out(Rcpp::NumericVector geno, 
		      int mrk, 
		      Rcpp::IntegerVector segreg_type, 
		      int n_ind, int verbose)
{
  int n_mar=((int)geno.size()/n_ind), k1, k2;
  int chunk=((n_mar*n_mar)-n_mar)/20, ct1=0, ct2=1, a1, a2, a3;
  NumericMatrix n(5,5);
  NumericVector r(8);
  NumericVector d1d2 =  NumericVector::create(0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00);
  NumericMatrix r1(n_mar,n_mar);
  NumericMatrix r2(n_mar,n_mar);
  NumericMatrix r3(n_mar,n_mar);
  NumericMatrix r4(n_mar,n_mar);
  NumericMatrix rm(4,n_mar);
  NumericMatrix lm(4,n_mar);
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
      R_CheckUserInterrupt(); /* check for ^C */
      if(mrk < 0) a3=(i+1);
      else a3=0;
      for(int j=a3; j  < n_mar; j++)
	{
	  ct1++;
    // Safe version using vectors with proper size
    std::vector<int> k_sub(n_ind);
    std::vector<int> k1_sub(n_ind);
        
    // Safe copying with bounds checking
    for(int k = 0; k < n_ind; k++) {
      if((i*n_ind + k) < geno.size() && (j*n_ind + k) < geno.size()) {
        k_sub[k] = geno[i*n_ind + k];
        k1_sub[k] = geno[j*n_ind + k];
      } else {
        Rcpp::stop("Index out of bounds in marker data access");
      }
    }
        
    // Add size verification
    if(k_sub.size() != n_ind || k1_sub.size() != n_ind) {
      Rcpp::stop("Inconsistent vector sizes in computation");
    }
    
	  std::fill(n.begin(), n.end(), 0);
	  std::fill(r.begin(), r.end(), 0);
	  for(int k=0; k < n_ind; k++)
	    {
	      if(k_sub[k]==1)
		{
		  if (k1_sub[k]==1) n(1,1)++;
		  else if (k1_sub[k]==2) n(1,2)++;
		  else if (k1_sub[k]==3) n(1,3)++;
		  else if (k1_sub[k]==4) n(1,4)++;
		  else n(0,0)++;
		}
	      else if(k_sub[k]==2)
		{
		  if (k1_sub[k]==1) n(2,1)++;
		  else if (k1_sub[k]==2) n(2,2)++;
		  else if (k1_sub[k]==3) n(2,3)++;
		  else if (k1_sub[k]==4) n(2,4)++;
		  else n(0,0)++;
		}
	      else if(k_sub[k]==3)
		{
		  if (k1_sub[k]==1) n(3,1)++;
		  else if (k1_sub[k]==2) n(3,2)++;
		  else if (k1_sub[k]==3) n(3,3)++;
		  else if (k1_sub[k]==4) n(3,4)++;
		  else n(0,0)++;
		}
	      else if(k_sub[k]==4)
		{
		  if (k1_sub[k]==1) n(4,1)++;
		  else if (k1_sub[k]==2) n(4,2)++;
		  else if (k1_sub[k]==3) n(4,3)++;
		  else if (k1_sub[k]==4) n(4,4)++;
		  else n(0,0)++;
		}
	      else n(0,0)++;
	    }
	  k1=segreg_type(i); k2=segreg_type(j);
	  switch(k1){
	  case 1:
	    switch(k2){
	    case 1:
	      r=rf_A_A(n, n_ind, n(0,0)); 	      /*Markers A - A */
	      break;
	    case 2:
	      r=rf_A_B1(n, n_ind, n(0,0));	      /*Markers A - B1 */
	      break;
	    case 3:
	      r=rf_A_B2(n, n_ind, n(0,0));	      /*Markers A - B2 */
	      break;
	    case 4:
	      r=rf_A_B3(n, n_ind, n(0,0));	      /*Markers A - B3 */
	      break;
	    case 5:
	      r=rf_A_C(n, n_ind, n(0,0));	      /*Markers A - C */
	      break;
	    case 6:
	      r=rf_A_D1(n, n_ind, n(0,0));	      /*Markers A - D1 */
	      break;
	    case 7:
	      r=rf_A_D2(n, n_ind, n(0,0));	      /*Markers A - D2 */
	      break;
	    }
	    break;
	  case 2:
	    switch(k2){
	    case 1:
	      n=transpose_counts(n);
	      r=rf_A_B1(n,n_ind, n(0,0));	      /*Markers B1 - A*/
	      break;
	    case 2:
	      r=rf_B1_B1(n,n_ind, n(0,0));	      /*Markers B1 - B1*/
	      break;
	    case 3:
	      r=rf_B1_B2(n,n_ind, n(0,0));	      /*Markers B1 - B2*/
	      break;
	    case 4:
	      r=rf_B1_B3(n,n_ind, n(0,0));	      /*Markers B1 - B3*/
	      break;
	    case 5:
	      r=rf_B1_C(n,n_ind, n(0,0));	      /*Markers B1 - c*/
	      break;
	    case 6:
	      r=rf_B1_D1(n,n_ind, n(0,0));	      /*Markers B1 - D1*/
	      break;
	    case 7:
	      r=rf_B1_D2(n,n_ind, n(0,0));	      /*Markers B1 - D2*/
	      break;
	    }
	    break;
	  case 3:
	    switch(k2){
	    case 1:
	      n=transpose_counts(n);
	      r=rf_A_B2(n,n_ind, n(0,0));	      /*Markers B2 - A*/
	      break;
	    case 2:
	      n=transpose_counts(n);
	      r=rf_B1_B2(n,n_ind, n(0,0));	      /*Markers B2 - B1*/
	      break;
	    case 3:
	      r=rf_B2_B2(n,n_ind, n(0,0));	      /*Markers B2 - B2*/
	      break;
	    case 4:
	      r=rf_B2_B3(n,n_ind, n(0,0));	      /*Markers B2 - B3*/
	      break;
	    case 5:
	      r=rf_B2_C(n,n_ind, n(0,0));	      /*Markers B2 - C*/
	      break;
	    case 6:
	      r=rf_B2_D1(n,n_ind, n(0,0));	      /*Markers B2 - D1*/
	      break;
	    case 7:
	      r=rf_B2_D2(n,n_ind, n(0,0));	      /*Markers B2 - D2*/
	      break;
	    }
	    break;
	  case 4:
	    switch(k2){
	    case 1:
	      n=transpose_counts(n);
	      r=rf_A_B3(n,n_ind, n(0,0));	      /*Markers B3 - A*/
	      break;
	    case 2:
	      n=transpose_counts(n);
	      r=rf_B1_B3(n,n_ind, n(0,0));	      /*Markers B3 - B1*/
	      break;
	    case 3:
	      n=transpose_counts(n);
	      r=rf_B2_B3(n,n_ind, n(0,0));	      /*Markers B3 - B2*/
	      break;
	    case 4:
	      r=rf_B3_B3(n,n_ind, n(0,0));	      /*Markers B3 - B3*/
	      break;
	    case 5:
	      r=rf_B3_C(n,n_ind, n(0,0));	      /*Markers B3 - C*/
	      break;
	    case 6:
	      r=rf_B3_D1(n,n_ind, n(0,0));	      /*Markers B3 - D1*/
	      break;
	    case 7:
	      r=rf_B3_D2(n,n_ind, n(0,0));	      /*Markers B3 - D2*/
	      break;
	    }
	    break;
	  case 5:
	    switch(k2){
	    case 1:
	      n=transpose_counts(n);
	      r=rf_A_C(n,n_ind, n(0,0));	      /*Markers C - A*/
	      break;
	    case 2:
	      n=transpose_counts(n);
	      r=rf_B1_C(n,n_ind, n(0,0));	      /*Markers C - B1*/
	      break;
	    case 3:
	      n=transpose_counts(n);
	      r=rf_B2_C(n,n_ind, n(0,0));	      /*Markers C - B2*/
	      break;
	    case 4:
	      n=transpose_counts(n);
	      r=rf_B3_C(n,n_ind, n(0,0));	      /*Markers C - B3*/
	      break;
	    case 5:
	      r=rf_C_C(n,n_ind, n(0,0));	      /*Markers C - C*/
	      break;
	    case 6:
	      r=rf_C_D1(n,n_ind, n(0,0));	      /*Markers C - D1*/
	      break;
	    case 7:
	      r=rf_C_D2(n,n_ind, n(0,0));	      /*Markers C - D2*/
	      break;
	    }
	    break;
	  case 6:
	    switch(k2){
	    case 1:
	      n=transpose_counts(n);
	      r=rf_A_D1(n,n_ind, n(0,0));	      /*Markers D1 - A*/
	      break;
	    case 2:
	      n=transpose_counts(n);
	      r=rf_B1_D1(n,n_ind, n(0,0));	      /*Markers D1 - B1*/
	      break;
	    case 3:
	      n=transpose_counts(n);
	      r=rf_B2_D1(n,n_ind, n(0,0));	      /*Markers D1 - B1*/
	      break;
	    case 4:
	      n=transpose_counts(n);
	      r=rf_B3_D1(n,n_ind, n(0,0));	      /*Markers D1 - B3*/
	      break;
	    case 5:
	      n=transpose_counts(n);
	      r=rf_C_D1(n,n_ind, n(0,0));	      /*Markers D1 - C*/
	      break;
	    case 6:
	      r=rf_D1_D1(n,n_ind, n(0,0));	      /*Markers D1 - D1*/
	      break;
	    case 7:
	      r= d1d2; //rep( NumericVector::get_na(), 8 );   /*Markers D1 - D2 - Impossible to compute*/
	      break;
	    }
	    break;
	  case 7:
	    switch(k2){
	    case 1:
	      n=transpose_counts(n);
	      r=rf_A_D2(n,n_ind, n(0,0));	      /*Markers D2 - A*/
	      break;
	    case 2:
	      n=transpose_counts(n);
	      r=rf_B1_D2(n,n_ind, n(0,0));	      /*Markers D2 -  B1*/
	      break;
	    case 3:
	      n=transpose_counts(n);
	      r=rf_B2_D2(n,n_ind, n(0,0));	      /*Markers D2 -  B2*/
	      break;
	    case 4:
	      n=transpose_counts(n);
	      r=rf_B3_D2(n,n_ind, n(0,0));	      /*Markers D2 -  B3*/
	      break;
	    case 5:
	      n=transpose_counts(n);
	      r=rf_C_D2(n,n_ind, n(0,0));	      /*Markers D2 -  C*/
	      break;
	    case 6:
	      r= d1d2;//rep( NumericVector::get_na(), 8 );   /*Markers D2 -  D1 - Impossible to compute*/
	      break;
	    case 7:
	      n=transpose_counts(n);
	      r=rf_D2_D2(n,n_ind, n(0,0));	      /*Markers D2 -  D2*/
	      break;
	    }
	    break;
	  }
	  for(int k=0; k < 4; k++)
	    {
	      if(r[k] < rf_TOL_min) r[k]=rf_TOL_min;
	      if(r[k] > rf_TOL_max) r[k]=rf_TOL_max;
	    }
	  if(mrk < 0)
	    {
	      r1(j,i)=r[0];
	      r2(j,i)=r[1];
	      r3(j,i)=r[2];
	      r4(j,i)=r[3];
	      r1(i,j)=r[4];
	      r2(i,j)=r[5];
	      r3(i,j)=r[6];
	      r4(i,j)=r[7];
	    }
	  else
	    {
	      rm(0,j)=r[0];
	      rm(1,j)=r[1];
	      rm(2,j)=r[2];
	      rm(3,j)=r[3];
	      lm(0,j)=r[4];
	      lm(1,j)=r[5];
	      lm(2,j)=r[6];
	      lm(3,j)=r[7];
	    }
	}
    }
  if(mrk < 0) 
    {
      if(verbose==1 && n_mar > 100) Rcpp::Rcout << "\t100%\n";
      List z = List::create(wrap(r1), wrap(r2), wrap(r3), wrap(r4));
      return(z);
    }
  else
    {
      List z = List::create(wrap(rm), wrap(lm));
      return(z);
    }
}
