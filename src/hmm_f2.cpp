/* Written by Marcelo Mollinari
   Adapted from hmm_main.c, hmm_f2.c and util.c (found in the R package qtl)
   copyright (c) 2001-10, Karl W Broman    

   These codes are under the GNU General Public License, version 3
   A copy of the GPL 3 is available at http://www.r-project.org/Licenses/GPL-3
*/


#include <Rcpp.h>
#include "hmm_f2.h"
using namespace Rcpp;
using namespace std;
#define THRESH 200.0

RcppExport SEXP est_hmm_f2(SEXP geno_R, SEXP rf_R, SEXP verbose_R){
  Rcpp::NumericMatrix geno = Rcpp::as<Rcpp::NumericMatrix>(geno_R);
  Rcpp::NumericVector rf = Rcpp::as<Rcpp::NumericVector>(rf_R);
  int verbose = Rcpp::as<int>(verbose_R);
  int n_mar = geno.nrow();
  int n_ind = geno.ncol();
  int n_gen = 4;
  int it, i, v, v2, j, j2, flag=0, maxit=1000;
  double error_prob = 0.00001;
  double tol=1e-6, s=0.0; 
  double loglik, curloglik; 
  NumericMatrix alpha(n_gen, n_mar);
  NumericMatrix beta(n_gen, n_mar);
  NumericMatrix gamma(n_gen, n_gen);
  NumericVector cur_rf(n_mar-1);
  NumericVector initf(4,0.25);

  NumericMatrix tr(n_gen, (n_mar-1)*n_gen);

  NumericMatrix em(6,4);
  em(0,0)=em(0,1)=em(0,2)=em(0,3)=1.0;
  em(1,0)=1.0-error_prob;
  em(1,1)=em(1,2)=em(1,3)=error_prob/2.0;
  em(2,1)=em(2,2)=1.0-error_prob;
  em(2,0)=em(2,0)=error_prob/2.0;
  em(3,3)=1.0-error_prob;
  em(3,0)=em(3,1)=em(3,2)=error_prob/2.0;
  em(4,0)=em(4,1)=em(4,2)=1.0-error_prob/2;
  em(4,3)=error_prob;
  em(5,1)=em(5,2)=em(5,3)=1.0-error_prob/2;;
  em(5,0)=error_prob;

  for(it=0; it<maxit; it++) {
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
	alpha(v,0) = initf(v)  * em(geno(0,i), v);
	beta(v,n_mar-1) = 1.0;
      }   
      /* forward-backward equations */
      for(j=1,j2=n_mar-2; j<n_mar; j++, j2--) {
	for(v=0; v<n_gen; v++) {
	  alpha(v,j) = alpha(0,j-1) * tr(0, (j-1)*n_gen+v);
	  beta(v,j2) = beta(0,j2+1) * tr(v, j2*4) * em(geno(j2+1,i),0);
	  for(v2=1; v2<n_gen; v2++) {
	    alpha(v,j) = alpha(v,j) + alpha(v2,j-1) * tr(v2,(j-1)*n_gen+v);
	    beta(v,j2) = beta(v,j2) + beta(v2,j2+1) * tr(v, j2*n_gen+v2)  * em(geno(j2+1,i),v2);
	  }
	  alpha(v,j) *= em(geno(j,i),v);
	}
      }
      for(j=0; j<n_mar-1; j++) {
	/* calculate gamma = log Pr(v1, v2, O) */
	for(v=0, s=0.0; v<n_gen; v++) {
	  for(v2=0; v2<n_gen; v2++) {
	    gamma(v,v2) = alpha(v,j) * beta(v2,j+1) * em(geno(j+1,i), v2)* tr(v, j*n_gen+v2);
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
      alpha(v,0) = initf(v) * em(geno(0,i), v);
    }
    // forward equations 
    for(j=1; j<n_mar; j++) {
      for(v=0; v<n_gen; v++) {
	alpha(v,j) = alpha(0,j-1) *
	  stepf_f2(1, v+1, rf(j-1));

	for(v2=1; v2<n_gen; v2++)
	  alpha(v,j) = alpha(v,j) + alpha(v2,j-1) *
	    stepf_f2(v2+1,v+1,rf(j-1));
	alpha(v,j) *= em(geno(j,i),v);
      }
    }
    curloglik = alpha(0,n_mar-1);
    for(v=1; v<n_gen; v++)
      curloglik = curloglik + alpha(v,n_mar-1);
    loglik += log(curloglik);
  }
  
  List z = List::create(wrap(rf), wrap(loglik));
  return(z);
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
