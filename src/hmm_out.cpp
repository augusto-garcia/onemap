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
 File: hmm_out.cpp
 Description: Set of functions to be used with software R
 Implements the methodology of Hidden Markov Models (HMM)
 to construct multipoint linkage maps in outcrossing species
 
 Written by Marcelo Mollinari
 Adapted from hmm_main.c, hmm_f2.c and util.c (found in the R package qtl)
 copyright (c) 2001-10, Karl W Broman
 
 Escola Superior de Agricultura "Luiz de Queiroz"
 Departamento de Genética - São Paulo, Brazil
 Contact: mmollina@usp.br
 First version in C: 03/02/2009
 First version in C++: 01/27/2016 
 Last update: 01/27/2016 
 */
#include <Rcpp.h>
#include "hmm_out.h"
using namespace Rcpp;
using namespace std;

//' Run HMM chains
//' 
//' @param geno_R genotypes
//' @param error_R genotypes probabilities
//' @param type_R marker types
//' @param phase_R phase estimated by 2-pts
//' @param rf_R recombination fraction
//' @param verbose_R logical to display or not the procedure
//' @param tol_R EM algorithm tolerance
//' 
//' @export
//' 
// [[Rcpp::export]]
RcppExport SEXP est_hmm_out(SEXP geno_R, SEXP error_R, SEXP type_R, SEXP phase_R, SEXP rf_R, SEXP verbose_R, SEXP tol_R){
  Rcpp::NumericMatrix geno = Rcpp::as<Rcpp::NumericMatrix>(geno_R);
  Rcpp::NumericMatrix error = Rcpp::as<Rcpp::NumericMatrix>(error_R);
  Rcpp::NumericVector type = Rcpp::as<Rcpp::NumericVector>(type_R);
  Rcpp::NumericVector phase = Rcpp::as<Rcpp::NumericVector>(phase_R);
  Rcpp::NumericVector rf = Rcpp::as<Rcpp::NumericVector>(rf_R);
  int verbose = Rcpp::as<int>(verbose_R);
  double tol = Rcpp::as<double>(tol_R);
  int n_mar = geno.nrow();
  int n_ind = geno.ncol();
  int n_gen = 4;
  int it, i, v, v2, j, j2, flag=0, maxit=1000;
  double s=0.0;
  double loglik, curloglik;
  NumericMatrix alpha(n_gen, n_mar);
  NumericMatrix beta(n_gen, n_mar);
  NumericMatrix gamma(n_gen, n_gen);
  NumericVector cur_rf(n_mar-1);
  NumericVector initf(4,0.25);
  NumericMatrix tr(n_gen, (n_mar-1)*n_gen);
  NumericMatrix probs(n_gen, (n_mar)*(n_ind));
  
  if(verbose) {
    //print initial estimates
    Rprintf("      ");
    for(int j=0; j<n_mar-1; j++) Rprintf("%.3lf ", rf[j]);
    Rprintf("\n");
  }
  // begin EM algorithm  
  for(it=0; it<maxit; it++) {
    for(j=0; j<n_mar-1; j++)
      cur_rf[j] = rf[j];
    //  1. If one wants to test similar orders, it is not necessary to
    //  allocate space for all markers again This is quite complicate to
    //  program, but could save a lot of processing time
    //  2. Consider genotyping errors for each marker in a similar way
    //  as was done for the trasition matrix
    for(j=0; j < ((n_mar-1)*n_gen); j++){
      for(i=0; i<n_gen; i++){
        tr(i,j)= step_out(i+1, (j%4)+1,  phase((int)j/4), cur_rf((int)j/4));
      }
    }
    std::fill(rf.begin(), rf.end(), 0);
    for(i=0; i<n_ind; i++) { // i = individual
      R_CheckUserInterrupt(); // check for ^C
      // initialize alpha and beta
      for(v=0; v<n_gen; v++) {
        alpha(v,0) = initf(v) * error(i*n_mar,v);
        beta(v,n_mar-1) = 1.0;
      }
      // forward-backward equations
      for(j=1,j2=n_mar-2; j<n_mar; j++, j2--) 
      {    
        for(v=0; v<n_gen; v++) {
          alpha(v,j) = alpha(0,j-1) * 
            tr(0, (j-1)*n_gen+v);     
          beta(v,j2) = beta(0,j2+1) * 
            tr(v,j2*n_gen)  * 
            error((j2+1)+(i*n_mar),0);
          for(v2=1; v2<n_gen; v2++) {
            alpha(v,j) = alpha(v,j) + alpha(v2,j-1) * 
              tr(v2,(j-1)*n_gen+v);
            beta(v,j2) = beta(v,j2) + beta(v2,j2+1) * 
              tr(v, j2*n_gen+v2) * 
              error((j2+1)+(i*n_mar), v2);
          }
          alpha(v,j) *= error(j + (i*n_mar),v);
        } 
      } 
      for(j=0; j<n_mar-1; j++) {
        // calculate gamma = log Pr(v1, v2, O)
        for(v=0, s=0.0; v<n_gen; v++) {
          for(v2=0; v2<n_gen; v2++) {
            gamma(v,v2) = alpha(v,j) * 
              beta(v2,j+1) * 
              error((j+1)+(i*n_mar),v2) * 
              tr(v, j*n_gen+v2);  
            if(v==0 && v2==0) s = gamma(v,v2);
            else s += gamma(v,v2);
          }
        }
        if(s==0) s = 0.000001; // avoiding NaN
        for(v=0; v<n_gen; v++) {
          for(v2=0; v2<n_gen; v2++) {
            rf(j) += nrec_out(v+1, v2+1, phase(j)) * gamma(v,v2)/s;
            if(ISNAN(rf(j))) {
              stop("Can not be estimated\n");
            }
          }
        }
      }
      
      /* Store genotypes probabilities*/
      for(j=0; j<n_mar; j++) {
        long double w = 0.0;
        for(v=0; v<n_gen; v++){
          w += alpha(v,j) * beta(v,j);
        }
        for(v=0; v<n_gen; v++){
          probs(v,j+i*n_mar) = (alpha(v,j) * beta(v,j)/w);
        }
      }
    } // loop over individuals
    
    // rescale
    for(j=0; j<n_mar-1; j++) {
      rf(j) /= (double)n_ind;
      if(rf(j) < tol/100.0) rf(j) = tol/100.0;
      else if(rf(j) > 0.5-tol/100.0) rf(j) = 0.5-tol/100.0;
    }
    // check convergence
    for(j=0, flag=0; j<n_mar-1; j++) {
      if(fabs(rf(j) - cur_rf(j)) > tol*(cur_rf(j)+tol*100.0)) {
        flag = 1;
        break;
      }
    }
    if(!flag) break;
  } // end EM algorithm
  if(flag) warning("Didn't converge!\n");
  
  // calculate log likelihood
  loglik = 0.0;
  for(i=0; i<n_ind; i++) { // i = individual
    // initialize alpha
    for(v=0; v<n_gen; v++)
      alpha(v,0) = initf(v) * error(i*n_mar,v);
    // forward equations
    for(j=1; j<n_mar; j++) {
      for(v=0; v<n_gen; v++) {
        alpha(v,j) = alpha(0,j-1) * 
          step_out(1, v+1, phase(j-1), rf(j-1));
        
        for(v2=1; v2<n_gen; v2++)
          alpha(v,j) = alpha(v,j) + alpha(v2,j-1) * 
            step_out(v2+1,v+1,phase(j-1),rf(j-1));
        alpha(v,j) *= error(j+(i*n_mar),v);
      }
    }
    curloglik = alpha(0,n_mar-1);
    for(v=1; v<n_gen; v++)
      curloglik = curloglik + alpha(v,n_mar-1);
    loglik += log(curloglik);
  }
  if(verbose) {
    //print final estimates
    Rprintf(" Number of iterations to converge: %4d \n", it+1);
    for(j=0; j<n_mar-1; j++) Rprintf(" %.3lf", rf(j));
    Rprintf("\n");
    
    Rprintf(" loglike: %10.4lf\n\n", loglik);
  }
  List z = List::create(wrap(rf), wrap(loglik), wrap(probs));
  return(z);
}

//Transition probability funtion for outcrossing species
double step_out(int gen1, int gen2, int phase, double rf)
{
  if(phase==1){/*CC*/
      switch(gen1) {
      case 1:
        switch(gen2) {
        case 1: return((1.0-rf)*(1.0-rf));
        case 2: return((1.0-rf)*rf);
        case 3: return((1.0-rf)*rf);
        case 4: return(rf*rf);
        }
      case 2:
        switch(gen2) {
        case 1: return((1.0-rf)*rf);
        case 2: return((1.0-rf)*(1.0-rf));
        case 3: return(rf*rf);
        case 4: return((1.0-rf)*rf);
        }
      case 3:
        switch(gen2) {
        case 1: return((1.0-rf)*rf);
        case 2: return(rf*rf);
        case 3: return((1.0-rf)*(1.0-rf));
        case 4: return((1.0-rf)*rf);
        }
      case 4:
        switch(gen2) {
        case 1: return(rf*rf);
        case 2: return((1.0-rf)*rf);
        case 3: return((1.0-rf)*rf);
        case 4: return((1.0-rf)*(1.0-rf));
        }
      }
    return(log(-1.0)); /* shouldn't get here */
  }
  else if(phase==2){/*CR*/
      switch(gen1) {
      case 1:
        switch(gen2) {
        case 1: return((1.0-rf)*rf);
        case 2: return((1.0-rf)*(1.0-rf));
        case 3: return(rf*rf);
        case 4: return((1.0-rf)*rf);
        }
      case 2:
        switch(gen2) {
        case 1: return((1.0-rf)*(1.0-rf));
        case 2: return((1.0-rf)*rf);
        case 3: return((1.0-rf)*rf);
        case 4: return(rf*rf);
        }
      case 3:
        switch(gen2) {
        case 1: return(rf*rf);
        case 2: return((1.0-rf)*rf);
        case 3: return((1.0-rf)*rf);
        case 4: return((1.0-rf)*(1.0-rf));
        }
      case 4:
        switch(gen2) {
        case 1: return((1.0-rf)*rf);
        case 2: return(rf*rf);
        case 3: return((1.0-rf)*(1.0-rf));
        case 4: return((1.0-rf)*rf);
        }
      }
    return(log(-1.0)); /* shouldn't get here */
  }
  else if(phase==3){
    switch(gen1) {
    case 1:
      switch(gen2) {
      case 1: return((1.0-rf)*rf);
      case 2: return(rf*rf);
      case 3: return((1.0-rf)*(1.0-rf));
      case 4: return((1.0-rf)*rf);
      }
    case 2:
      switch(gen2) {
      case 1: return(rf*rf);
      case 2: return((1.0-rf)*rf);
      case 3: return((1.0-rf)*rf);
      case 4: return((1.0-rf)*(1.0-rf));
      }
    case 3:
      switch(gen2) {
      case 1: return((1.0-rf)*(1.0-rf));
      case 2: return((1.0-rf)*rf);
      case 3: return((1.0-rf)*rf);
      case 4: return(rf*rf);
      }
    case 4:
      switch(gen2) {
      case 1: return((1.0-rf)*rf);
      case 2: return((1.0-rf)*(1.0-rf));
      case 3: return(rf*rf);
      case 4: return((1.0-rf)*rf);
      }
    }
    return(log(-1.0)); /* shouldn't get here */
  }
  else if(phase==4){
    switch(gen1) {
    case 1:
      switch(gen2) {
      case 1: return(rf*rf);
      case 2: return((1.0-rf)*rf);
      case 3: return((1.0-rf)*rf);
      case 4: return((1.0-rf)*(1.0-rf));
      }
    case 2:
      switch(gen2) {
      case 1: return((1.0-rf)*rf);
      case 2: return(rf*rf);
      case 3: return((1.0-rf)*(1.0-rf));
      case 4: return((1.0-rf)*rf);
      }
    case 3:
      switch(gen2) {
      case 1: return((1.0-rf)*rf);
      case 2: return((1.0-rf)*(1.0-rf));
      case 3: return(rf*rf);
      case 4: return((1.0-rf)*rf);
      }
    case 4:
      switch(gen2) {
      case 1: return((1.0-rf)*(1.0-rf));
      case 2: return((1.0-rf)*rf);
      case 3: return((1.0-rf)*rf);
      case 4: return(rf*rf);
      }
    }
    return(0.0); /* shouldn't get here */
  }
  else return(0.0); /* shouldn't get here */
}

//Expected recombination events for outcrossing species
double nrec_out(int gen1, int gen2, int phase)
{
  if(phase==1){
    switch(gen1) {
    case 1:
      switch(gen2) {
      case 1: return(0);
      case 2: return(0.5);
      case 3: return(0.5);
      case 4: return(1);
      }
    case 2:
      switch(gen2) {
      case 1: return(0.5);
      case 2: return(0);
      case 3: return(1);
      case 4: return(0.5);
      }
    case 3:
      switch(gen2) {
      case 1: return(0.5);
      case 2: return(1);
      case 3: return(0);
      case 4: return(0.5);
      }
    case 4:
      switch(gen2) {
      case 1: return(1);
      case 2: return(0.5);
      case 3: return(0.5);
      case 4: return(0);
      }
    }
    return(log(-1.0)); /* shouldn't get here */
  }
  else if(phase==2){
    switch(gen1) {
    case 1:
      switch(gen2) {
      case 1: return(0.5);
      case 2: return(0);
      case 3: return(1);
      case 4: return(0.5);
      }
    case 2:
      switch(gen2) {
      case 1: return(0);
      case 2: return(0.5);
      case 3: return(0.5);
      case 4: return(1);
      }
    case 3:
      switch(gen2) {
      case 1: return(1);
      case 2: return(0.5);
      case 3: return(0.5);
      case 4: return(0);
      }
    case 4:
      switch(gen2) {
      case 1: return(0.5);
      case 2: return(1);
      case 3: return(0);
      case 4: return(0.5);
      }
    }
    return(log(-1.0)); /* shouldn't get here */
  }
  else if(phase==3){
    switch(gen1) {
    case 1:
      switch(gen2) {
      case 1: return(0.5);
      case 2: return(1);
      case 3: return(0);
      case 4: return(0.5);
      }
    case 2:
      switch(gen2) {
      case 1: return(1);
      case 2: return(0.5);
      case 3: return(0.5);
      case 4: return(0);
      }
    case 3:
      switch(gen2) {
      case 1: return(0);
      case 2: return(0.5);
      case 3: return(0.5);
      case 4: return(1);
      }
    case 4:
      switch(gen2) {
      case 1: return(0.5);
      case 2: return(0);
      case 3: return(1);
      case 4: return(0.5);
      }
    }
    return(log(-1.0)); /* shouldn't get here */
  }
  else if(phase==4){
    switch(gen1) {
    case 1:
      switch(gen2) {
      case 1: return(1);
      case 2: return(0.5);
      case 3: return(0.5);
      case 4: return(0);
      }
    case 2:
      switch(gen2) {
      case 1: return(0.5);
      case 2: return(1);
      case 3: return(0);
      case 4: return(0.5);
      }
    case 3:
      switch(gen2) {
      case 1: return(0.5);
      case 2: return(0);
      case 3: return(1);
      case 4: return(0.5);
      }
    case 4:
      switch(gen2) {
      case 1: return(0);
      case 2: return(0.5);
      case 3: return(0.5);
      case 4: return(1);
      }
    }
    return(log(-1.0)); /* shouldn't get here */
  }
  else return(log(-1.0)); /* shouldn't get here */
}

