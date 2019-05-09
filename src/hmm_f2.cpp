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
  File: hmm_f2.cpp
  Description: Set of functions to be used with software R
               Implements the methodology of Hidden Markov Models (HMM)
	       to construct multipoint linkage maps in f2 crosses

  Written by Marcelo Mollinari
  Adapted from hmm_main.c, hmm_f2.c and util.c (found in the R package qtl)
  copyright (c) 2001-10, Karl W Broman                                

  Escola Superior de Agricultura "Luiz de Queiroz"
  Departamento de Genética - São Paulo, Brazil
  Contact: mmollina@usp.br
  First version: 10/2015
  Last update: 12/2015
*/

#include <Rcpp.h>
#include "hmm_f2.h"
using namespace Rcpp;
using namespace std;
#define THRESH 200.0


/**********************************************************************
 * 
 * est_hmm_f2
 *
 * This function re-estimates the genetic map for a chromosome in a f2 cross
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

RcppExport SEXP est_hmm_f2(SEXP geno_R, SEXP error_R, SEXP rf_R, SEXP verbose_R, SEXP tol_R){
  Rcpp::NumericMatrix geno = Rcpp::as<Rcpp::NumericMatrix>(geno_R);
  Rcpp::NumericMatrix error = Rcpp::as<Rcpp::NumericMatrix>(error_R);
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

  /*NumericMatrix em(6,4);
  em(0,0)=em(0,1)=em(0,2)=em(0,3)=1.0;
  em(1,0)=1.0-error_prob;
  em(1,1)=em(1,2)=em(1,3)=error_prob/3.0;
  em(2,1)=em(2,2)=1.0-error_prob;
  em(2,0)=em(2,3)=error_prob/2.0;
  em(3,3)=1.0-error_prob;
  em(3,0)=em(3,1)=em(3,2)=error_prob/3.0;
  em(4,0)=em(4,1)=em(4,2)=1.0-error_prob/3.0;
  em(4,3)=error_prob;
  em(5,1)=em(5,2)=em(5,3)=1.0-error_prob/3.0;;
  em(5,0)=error_prob;*/
  if(verbose) {
    /* print initial estimates */
    Rprintf("      "); 
    for(int j=0; j<n_mar-1; j++) Rprintf("%.3lf ", rf[j]);
    Rprintf("\n"); 
  }
  for(it=0; it<maxit; it++) {
    /* if(verbose) {
      Rprintf("\n"); 
      Rprintf("it: %4d ", it+1);
      for(j=0; j<n_mar-1; j++) Rprintf("%.3lf ", rf[j]);
    }
    */
    //Rcpp::Rcout << "it: " << it << "\n";
    for(j=0; j<n_mar-1; j++) {
      cur_rf[j] = rf[j];
      rf[j] = 0.0;
    }
    /*
      1. If one wants to test similar orders, it is not necessary to
      allocate space for all markers again This is quite complicate to
      program, but could save a lot of processing time 

      2. Consider genotyping errors for each marker in a similar way
      as was done for the trasition matrix
    */
    for(j=0; j < ((n_mar-1)*n_gen); j++){
      for(i=0; i<n_gen; i++){
	tr(i,j)= stepf_f2(i+1, (j%4)+1, cur_rf(j/4));
      }
    }

    std::fill(rf.begin(), rf.end(), 0);
    for(i=0; i<n_ind; i++) { /* i = individual */
      R_CheckUserInterrupt(); /* check for ^C */
      /* initialize alpha and beta */
      for(v=0; v<n_gen; v++) {
	alpha(v,0) = initf(v)  * emit_inb(geno(0,i), v+1, error(0,i));
	beta(v,n_mar-1) = 1.0;
      }   
      /* forward-backward equations */
      for(j=1,j2=n_mar-2; j<n_mar; j++, j2--) {
	for(v=0; v<n_gen; v++) {
	  alpha(v,j) = alpha(0,j-1) * tr(0, (j-1)*n_gen+v);
	  beta(v,j2) = beta(0,j2+1) * tr(v, j2*4) * emit_inb(geno(j2+1,i),1, error(j2+1,i));
	  for(v2=1; v2<n_gen; v2++) {
	    alpha(v,j) = alpha(v,j) + alpha(v2,j-1) * tr(v2,(j-1)*n_gen+v);
	    beta(v,j2) = beta(v,j2) + beta(v2,j2+1) * tr(v, j2*n_gen+v2)  * emit_inb(geno(j2+1,i),v2+1,error(j2+1,i));
	  }
	  alpha(v,j) *= emit_inb(geno(j,i),v+1, error(j,i));
	}
      }
      for(j=0; j<n_mar-1; j++) {
	/* calculate gamma = log Pr(v1, v2, O) */
	for(v=0, s=0.0; v<n_gen; v++) {
	  for(v2=0; v2<n_gen; v2++) {
	    gamma(v,v2) = alpha(v,j) * beta(v2,j+1) * emit_inb(geno(j+1,i), v2+1, error(j+1,i))* tr(v, j*n_gen+v2);
	    if(v==0 && v2==0) s = gamma(v,v2);
	    else s = s + gamma(v,v2);
	  }
	}
	for(v=0; v<n_gen; v++) {
	  for(v2=0; v2<n_gen; v2++) {
	    rf(j) += nrecf_f2(v+1,v2+1) * gamma(v,v2)/s;
	  }
	}
      }
    } /* loop over individuals */
    /* rescale */
    for(j=0; j<n_mar-1; j++) {
      rf(j) /= (double)n_ind;
      if(rf(j) < tol/1000.0) rf(j) = tol/1000.0;
      else if(rf(j) > 0.5-tol/1000.0) rf(j) = 0.5-tol/1000.0;
    }
    /* check convergence */    
    for(j=0, flag=0; j<n_mar-1; j++) {
      if(fabs(rf(j) - cur_rf(j)) > tol*(cur_rf(j)+tol*100.0)) {
	flag = 1;
	break;
      }
    }
    if(!flag) break;

  } /* end EM algorithm */

  //if(flag) warning("Didn't converge!\n");

  // calculate log likelihood
  loglik = 0.0;
  
  for(i=0; i<n_ind; i++) { // i = individual 
    // initialize alpha 
    for(v=0; v<n_gen; v++) {
      alpha(v,0) = initf(v) * emit_inb(geno(0,i), v+1, error(0,i));
    }
    // forward equations 
    for(j=1; j<n_mar; j++) {
      for(v=0; v<n_gen; v++) {
	alpha(v,j) = alpha(0,j-1) *
	  stepf_f2(1, v+1, rf(j-1));

	for(v2=1; v2<n_gen; v2++)
	  alpha(v,j) = alpha(v,j) + alpha(v2,j-1) *
	    stepf_f2(v2+1,v+1,rf(j-1));
	alpha(v,j) *= emit_inb(geno(j,i),v+1, error(j,i));
      }
    }
    curloglik = alpha(0,n_mar-1);
    for(v=1; v<n_gen; v++)
      curloglik = curloglik + alpha(v,n_mar-1);
    loglik += log(curloglik);
  }

  if(verbose) {
    /* print final estimates */
    Rprintf(" Number of iterations to converge: %4d \n", it+1);
    for(j=0; j<n_mar-1; j++) Rprintf(" %.3lf", rf[j]);
    Rprintf("\n");
    
    Rprintf(" loglike: %10.4lf\n\n", loglik);
  }  

  List z = List::create(wrap(rf), wrap(loglik));
  return(z);
}

//Emission function for inbred species

double emit_inb(int obs_gen, int true_gen, double error)
{
  /* Format: 6 x 4
     column(true gen):            A       AB    BA      B
     row (obs gen): missing(0)    1       1      1      1
     A(1)      1-e     e/3    e/3    e/3
     H(2)      e/2     1-e    1-e    e/2
     B(3)      e/3     e/3    e/3    1-e
     D(4)     1-e/3    1-e/3  1-e/3   e
     C(5)       e      1-e/3  1-e/3  1-e/3
  */
  switch(obs_gen){
  case 0: return(1.0);
  case 1:
    switch(true_gen){
    case 1: return(1.0-error);
    case 2: case 3: case 4: return(error/3.0);
    }
  case 2:
    switch(true_gen){
    case 2: case 3: return(1.0-error);
    case 1: case 4: return(error/2.0);
    }
  case 3:
    switch(true_gen){
    case 4: return(1.0-error);
    case 1: case 2: case 3: return(error/3.0);
    }
  case 4:
    switch(true_gen){
    case 1: case 2: case 3: return(1 - error/3.0);
    case 4: return(error);
    }
  case 5:
    switch(true_gen){
    case 1: return(error);
    case 2: case 3: case 4: return(1 - error/3.0);
    }    
  }
}


double stepf_f2(int gen1, int gen2, double rf)
{
  if(gen1==gen2) return((1.0-rf)*(1.0-rf));
  else if((gen1+gen2)==5) return(rf*rf);
  else return((1.0-rf)*rf);
}

double nrecf_f2(int gen1, int gen2)
{
  if(gen1==gen2) return(0.0);
  else if((gen1+gen2)==5) return(1.0);
  else return(0.5);
}
