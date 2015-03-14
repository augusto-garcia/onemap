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
  File: hmm_out.c
  Description: Set of functions to be used with software R
               Implements the methodology of Hidden Markov Models (HMM)
	       to construct multipoint linkage maps in outcrossing species

  Written by Marcelo Mollinari
  Adapted from hmm_main.c, hmm_f2.c and util.c (found in the R package qtl)
  copyright (c) 2001-10, Karl W Broman                                

  Escola Superior de Agricultura "Luiz de Queiroz"
  Departamento de Genética - São Paulo, Brazil
  Contact: mmollina@usp.br
  First version: 03/02/2009
  Last update: 03/02/2009
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h> 
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#define THRESH 200.0
#define TOL       1.0e-12
#define LN_05    -0.6931471805599453 /* natural log of 0.5 */
#define LN_025   -1.3862943611198906 /* natural log of 0.25 */
#define LN_0125  -2.0794415416798357 /* natural log of 0.125 */
#define LN_2      0.6931471805599453 /* natural log of 2 */
#define LN_3      1.09861228866811   /* natoral log of 3 */
                  
/**********************************************************************
 * 
 * addlog
 *
 * Calculate addlog(a,b) = log[exp(a) + exp(b)]
 *
 * This makes use of the function log1p(x) = log(1+x) provided
 * in R's math library.   
 *
 **********************************************************************/
double addlog(double a, double b)
{
  if(b > a + THRESH) return(b);
  else if(a > b + THRESH) return(a);
  else return(a + log1p(exp(b-a)));
}
		       

/**********************************************************************
 * 
 * reorg_geno
 *
 * Reorganize the genotype data so that it is a doubly indexed array
 * rather than a single long vector
 *
 * Afterwards, geno indexed like Geno[mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_geno(int n_ind, int n_pos, int *geno, int ***Geno)
{
  int i;

  *Geno = (int **)R_alloc(n_pos, sizeof(int *));

  (*Geno)[0] = geno;
  for(i=1; i< n_pos; i++) 
    (*Geno)[i] = (*Geno)[i-1] + n_ind;

}

/**********************************************************************
 * 
 * reorg_genoprob
 *
 * Reorganize the genotype probability data so that it is a triply 
 * indexed array rather than a single long vector
 *
 * Afterwards, genoprob indexed like Genoprob[gen][mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_genoprob(int n_ind, int n_pos, int n_gen, 
		    double *genoprob, double ****Genoprob)
{
  int i, j;
  double **a;

  *Genoprob = (double ***)R_alloc(n_gen, sizeof(double **));

  a = (double **)R_alloc(n_pos*n_gen, sizeof(double *));

  (*Genoprob)[0] = a;
  for(i=1; i< n_gen; i++) 
    (*Genoprob)[i] = (*Genoprob)[i-1]+n_pos;
  
  for(i=0; i<n_gen; i++) 
    for(j=0; j<n_pos; j++) 
      (*Genoprob)[i][j] = genoprob + i*n_ind*n_pos + j*n_ind;
}

/**********************************************************************
 * 
 * allocate_alpha
 *
 * Allocate space for alpha and beta matrices
 *
 * Afterwards, indexed like alpha[gen][mar]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_alpha(int n_pos, int n_gen, double ***alpha)
{
  int i;

  *alpha = (double **)R_alloc(n_gen, sizeof(double *));

  (*alpha)[0] = (double *)R_alloc(n_gen*n_pos, sizeof(double));

  for(i=1; i< n_gen; i++) 
    (*alpha)[i] = (*alpha)[i-1] + n_pos;
}


/**********************************************************************
 * 
 * allocate_double
 *
 * Allocate space for a vector of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_double(int n, double **vector)
{
  *vector = (double *)R_alloc(n, sizeof(double));
}

/**********************************************************************
 * 
 * allocate_int
 *
 * Allocate space for a vector of ints
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_int(int n, int **vector)
{
  *vector = (int *)R_alloc(n, sizeof(int));
}

/**********************************************************************
 * 
 * allocate_dmatrix
 *
 * Allocate space for a matrix of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_dmatrix(int n_row, int n_col, double ***matrix)
{
  int i;

  *matrix = (double **)R_alloc(n_row, sizeof(double *));
  
  (*matrix)[0] = (double *)R_alloc(n_col*n_row, sizeof(double));

  for(i=1; i<n_row; i++) 
    (*matrix)[i] = (*matrix)[i-1]+n_col;
}

/**********************************************************************
 * 
 * allocate_imatrix
 *
 * Allocate space for a matrix of ints
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_imatrix(int n_row, int n_col, int ***matrix)
{
  int i;

  *matrix = (int **)R_alloc(n_row, sizeof(int *));
  
  (*matrix)[0] = (int *)R_alloc(n_col*n_row, sizeof(int));

  for(i=1; i<n_row; i++) 
    (*matrix)[i] = (*matrix)[i-1]+n_col;
}


double init_outbred(int true_gen)
{
  return(LN_025);
}

double emit_outbred(int obs_gen, int true_gen, double error_prob, int mark_type)
{  
  /*Notation similar to Wu et. al.(2002)
    
  Marker type: A (1:1:1:1) (m=1 in the codes)
  Parental genotype codification (not used in the codes)
  
    A.1        A.2        A.3        A.4 (in the article)
  
  P1  P2     P1  P2     P1  P2     P1  P2      
  -a- -c-    -a- -a-    -a- -c-    -a- -b- 
     X    or    X    or    X    or    X   
  -b- -d-    -b- -c-    -b- -o-    -o- -o-
  
  offspring genotype codification (observed)
  missing: 0
  ac;  a; ac; ab: 1
  ad; ac;  a;  a: 2
  bc; ba; bc;  b: 3
  bd; bc;  b;  o: 4
  */ 
  switch(mark_type){
  case 1: /*A*/
    switch(obs_gen){
    case 0: return(0.0);
    case 1: case 2: case 3: case 4:
	if(obs_gen==true_gen) return(log(1-error_prob));
	else return(log(error_prob)-LN_3);
    }
    return(0.0);/* shouldn't get here */
    
  /*
    Marker type: B (2:1:1) (m=2, m=3, m=4 in the codes)
    Parental genotype codification (do not used in the codes)
    
      B.1        B.2        B.3     (in the article)
    
    P1  P2     P1  P2     P1  P2
    -a- -a-    -a- -a-    -a- -a-
       X    or    X    or    X   
    -b- -o-    -o- -b-    -b- -b-
    
    offspring genotype codification (observed)
    
    VERY IMPORTANT: The proportion here is 2:1:1, unlike 1:2:1, usualy used!!!!!  
    missing: 0
    a  a  ab : 1
    ab ab a  : 2 
    b  b  b  : 3
  */
    
  case 2: /*B.1*/ 
    switch(obs_gen){
    case 0: return(0.0);
    case 1: 
      switch(true_gen){
      case 1: case 2: return(log(1-error_prob));
      case 3: case 4: return(log(error_prob)-LN_2);
      }
    case 2:
      switch(true_gen){
      case 3: return(log(1-error_prob));
      case 1: case 2: case 4: return(log(error_prob)-LN_3);
      }
    case 3:
      switch(true_gen){
      case 4: return(log(1-error_prob));
      case 1: case 2: case 3: return(log(error_prob)-LN_3);
      }
    }
    return(0.0);/* shouldn't get here */
  
  case 3:/*B.2*/ 
    switch(obs_gen){
    case 0: return(0.0);
    case 1: 
      switch(true_gen){
      case 1: case 3: return(log(1-error_prob));
      case 2: case 4: return(log(error_prob)-LN_2);
      }
    case 2:
      switch(true_gen){
      case 2: return(log(1-error_prob));
      case 1: case 3: case 4: return(log(error_prob)-LN_3);
      }
    case 3:
      switch(true_gen){
      case 4: return(log(1-error_prob));
      case 1: case 2: case 3: return(log(error_prob)-LN_3);
      }
    }
    return(0.0);/* shouldn't get here */

  case 4: /*B.3*/ 
    switch(obs_gen){
    case 0: return(0.0);
    case 1: 
      switch(true_gen){
      case 1: return(log(1-error_prob));
      case 2: case 3: case 4: return(log(error_prob)-LN_3);
      }
    case 2:
      switch(true_gen){
      case 2: case 3: return(log(1-error_prob));
      case 1: case 4: return(log(error_prob)-LN_2);
     }
    case 3:
      switch(true_gen){
      case 4: return(log(1-error_prob));
      case 1: case 2: case 3: return(log(error_prob)-LN_3);
      }
    }
    return(0.0);/* shouldn't get here */ 
 

/*
    Marker type: C (3:1) (m=5, in the codes)
    Parental genotype codification (do not used in the codes)
    
      C   (in the article)
    
    P1  P2    
    -a- -a-
       X   
    -o- -o-
    
    offspring genotype codification (observed)
    missing: 0
    a: 1
    o: 2 
  */
    case 5:  /*C*/ 
      switch(obs_gen){
      case 0: return(0); 
      case 1: 
	if(true_gen==4) return(log(error_prob));
	else return(log(1-error_prob));
      case 2:
	if(true_gen==4) return(log(1-error_prob));
	else return(log(error_prob)-LN_3);
      }
      return(0.0);/* shouldn't get here */
  
/*
    Marker type: D (1:1) (m=6,m=7 in the codes)
    Parental genotype codification (do not used in the codes)
    
      D.1   (in the article)
    
    P1  P2    
    -a- -c-
       X   ...more 4 types, see on the article
    -b- -c-

      D.2   (in the article)
    
    P1  P2    
    -c- -a-
       X   ...more 4 types, see on the article
    -c- -b-

    offspring genotype codification (observed)
    missing: 0
    D.1
    ac: 1
    bc: 2 

    D.2
    ac: 2
    bc: 1 

  */
    case 6:  /*D.1*/ 
      switch(obs_gen){
      case 0: return(0);
      case 1: 
	if(true_gen==1||true_gen==2) return(log(1-error_prob));
	else return(log(error_prob)-LN_2);
      case 2:
	if(true_gen==3||true_gen==4) return(log(1-error_prob));
	else return(log(error_prob)-LN_2);
      }
      return(0.0);/* shouldn't get here */
  
    case 7:  /*D.2*/ 
      switch(obs_gen){
      case 0: return(0);
      case 1: 
	if(true_gen==1||true_gen==3) return(log(1-error_prob));
	else return(log(error_prob)-LN_2);
      case 2:
	if(true_gen==2||true_gen==4) return(log(1-error_prob));
	else return(log(error_prob)-LN_2);
      }
      return(0.0);/* shouldn't get here */

  }
  return(0.0);/* shouldn't get here */
}

double step_outbred(int gen1, int gen2, int phase, double rf)
{
  if(phase==1){/*CC*/
    switch(gen1) {
    case 1:
      switch(gen2) {
      case 1: return(2.0*log(1.0-rf));
      case 2: return(log(1.0-rf)+log(rf));
      case 3: return(log(1.0-rf)+log(rf));
      case 4: return(2.0*log(rf));
      }
    case 2:
      switch(gen2) {
      case 1: return(log(1.0-rf)+log(rf));
      case 2: return(2.0*log(1.0-rf));
      case 3: return(2.0*log(rf));
      case 4: return(log(1.0-rf)+log(rf));
      }
    case 3:
      switch(gen2) {
      case 1: return(log(1.0-rf)+log(rf));
      case 2: return(2.0*log(rf));
      case 3: return(2.0*log(1.0-rf));
      case 4: return(log(1.0-rf)+log(rf));
      }
    case 4:
      switch(gen2) {
      case 1: return(2.0*log(rf));
      case 2: return(log(1.0-rf)+log(rf));
      case 3: return(log(1.0-rf)+log(rf));
      case 4: return(2.0*log(1.0-rf));
      }
    }
    return(log(-1.0)); /* shouldn't get here */
  }
  else if(phase==2){/*CR*/
    switch(gen1) {
    case 1:
      switch(gen2) {
      case 1: return(log(1.0-rf)+log(rf));
      case 2: return(2.0*log(1.0-rf));
      case 3: return(2.0*log(rf));
      case 4: return(log(1.0-rf)+log(rf));
      }
    case 2:
      switch(gen2) {
      case 1: return(2.0*log(1.0-rf));
      case 2: return(log(1.0-rf)+log(rf));
      case 3: return(log(1.0-rf)+log(rf));
      case 4: return(2.0*log(rf));
      }
    case 3:
      switch(gen2) {
      case 1: return(2.0*log(rf));
      case 2: return(log(1.0-rf)+log(rf));
      case 3: return(log(1.0-rf)+log(rf));
      case 4: return(2.0*log(1.0-rf));
      }
    case 4:
      switch(gen2) {
      case 1: return(log(1.0-rf)+log(rf));
      case 2: return(2.0*log(rf));
      case 3: return(2.0*log(1.0-rf));
      case 4: return(log(1.0-rf)+log(rf));
      }
    }
    return(log(-1.0)); /* shouldn't get here */
  }
  else if(phase==3){
    switch(gen1) {
    case 1:
      switch(gen2) {
      case 1: return(log(1.0-rf)+log(rf));
      case 2: return(2.0*log(rf));
      case 3: return(2.0*log(1.0-rf));
      case 4: return(log(1.0-rf)+log(rf));
      }
    case 2:
      switch(gen2) {
      case 1: return(2.0*log(rf));
      case 2: return(log(1.0-rf)+log(rf));
      case 3: return(log(1.0-rf)+log(rf));
      case 4: return(2.0*log(1.0-rf));
      }
    case 3:
      switch(gen2) {
      case 1: return(2.0*log(1.0-rf));
      case 2: return(log(1.0-rf)+log(rf));
      case 3: return(log(1.0-rf)+log(rf));
      case 4: return(2.0*log(rf));
      }
    case 4:
      switch(gen2) {
      case 1: return(log(1.0-rf)+log(rf));
      case 2: return(2.0*log(1.0-rf));
      case 3: return(2.0*log(rf));
      case 4: return(log(1.0-rf)+log(rf));
      }
    }
    return(log(-1.0)); /* shouldn't get here */
  }
  else if(phase==4){
    switch(gen1) {
    case 1:
      switch(gen2) {
      case 1: return(2.0*log(rf));
      case 2: return(log(1.0-rf)+log(rf));
      case 3: return(log(1.0-rf)+log(rf));
      case 4: return(2.0*log(1.0-rf));
      }
    case 2:
      switch(gen2) {
      case 1: return(log(1.0-rf)+log(rf));
      case 2: return(2.0*log(rf));
      case 3: return(2.0*log(1.0-rf));
      case 4: return(log(1.0-rf)+log(rf));
      }
    case 3:
      switch(gen2) {
      case 1: return(log(1.0-rf)+log(rf));
      case 2: return(2.0*log(1.0-rf));
      case 3: return(2.0*log(rf));
      case 4: return(log(1.0-rf)+log(rf));
      }
    case 4:
      switch(gen2) {
      case 1: return(2.0*log(1.0-rf));
      case 2: return(log(1.0-rf)+log(rf));
      case 3: return(log(1.0-rf)+log(rf));
      case 4: return(2.0*log(rf));
      }
    }
    return(log(-1.0)); /* shouldn't get here */
  }
  else return(log(-1.0)); /* shouldn't get here */
}
  
  


double nrec_outbred(int gen1, int gen2, int phase){
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




/**********************************************************************
 * 
 * est_map
 *
 * This function re-estimates the genetic map for a chromosome
 *
 * n_ind        Number of individuals
 *
 * n_mar        Number of markers 
 *
 * n_gen        Number of different genotypes
 *
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 * error_prob   Genotyping error probability
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 * nrecf1       Function returning number of recombinations associated
 *              with (g_1, g_2)
 *
 * loglik       Loglik at final estimates of recombination fractions
 *
 * maxit        Maximum number of iterations to perform
 * 
 * tol          Tolerance for determining convergence
 * 
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void est_map(int n_ind, int n_mar, int *type, int *phase, int n_gen, int *geno, double *rf, 
	     double error_prob, double initf(int), 
	     double emitf(int, int, double, int),
	     double stepf(int, int, int, double), 
	     double nrec(int, int,int), 
	     double *loglik, int maxit, double tol, 
	     int verbose)
{
  int i, j, j2, v, v2, it, flag=0, **Geno;
  double s, **alpha, **beta, **gamma, *cur_rf;
  double curloglik;
  
  /* allocate space for beta and reorganize geno */
  reorg_geno(n_ind, n_mar, geno, &Geno);
  allocate_alpha(n_mar, n_gen, &alpha);
  allocate_alpha(n_mar, n_gen, &beta);
  allocate_dmatrix(n_gen, n_gen, &gamma);
  allocate_double(n_mar-1, &cur_rf);

  if(verbose) {
    /* print initial estimates */
    Rprintf("      "); 
    for(j=0; j<n_mar-1; j++) Rprintf("%.3lf ", rf[j]);
    Rprintf("\n"); 
  }

  /* begin EM algorithm */
  for(it=0; it<maxit; it++) {
    
    /*if(verbose) {
      Rprintf("\n"); 
      Rprintf("it: %4d ", it+1);
      for(j=0; j<n_mar-1; j++) Rprintf("%.3lf ", rf[j]);
      }
    */

    for(j=0; j<n_mar-1; j++) {
      cur_rf[j] = rf[j];
      rf[j] = 0.0;
    }

    for(i=0; i<n_ind; i++) { /* i = individual */

      /* initialize alpha and beta */
      for(v=0; v<n_gen; v++) {
	alpha[v][0] = initf(v+1) + emitf(Geno[0][i], v+1, error_prob, type[0]);
	beta[v][n_mar-1] = 0.0;
      }

      /* forward-backward equations */
      for(j=1,j2=n_mar-2; j<n_mar; j++, j2--) {
	
	for(v=0; v<n_gen; v++) {
	  alpha[v][j] = alpha[0][j-1] + stepf(1, v+1, phase[j-1], cur_rf[j-1]);
	  
	  beta[v][j2] = beta[0][j2+1] + stepf(v+1,1,phase[j2],cur_rf[j2]) + 
	    emitf(Geno[j2+1][i],1,error_prob,type[j2+1]);
	  
	  for(v2=1; v2<n_gen; v2++) {
	    alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
				 stepf(v2+1,v+1,phase[j-1],cur_rf[j-1]));
	    beta[v][j2] = addlog(beta[v][j2], beta[v2][j2+1] + 
				 stepf(v+1,v2+1,phase[j2],cur_rf[j2]) +
				 emitf(Geno[j2+1][i],v2+1,error_prob,type[j2+1]));
	  }
	  
	  alpha[v][j] += emitf(Geno[j][i],v+1,error_prob,type[j]);
		 
	}

      }

      for(j=0; j<n_mar-1; j++) {

	/* calculate gamma = log Pr(v1, v2, O) */
	for(v=0, s=0.0; v<n_gen; v++) {
	  for(v2=0; v2<n_gen; v2++) {
	    gamma[v][v2] = alpha[v][j] + beta[v2][j+1] + 
	      emitf(Geno[j+1][i], v2+1, error_prob, type[j+1]) +
	      stepf(v+1, v2+1, phase[j], cur_rf[j]);

	    if(v==0 && v2==0) s = gamma[v][v2];
	    else s = addlog(s, gamma[v][v2]);
	  }
	}

	for(v=0; v<n_gen; v++) {
	  for(v2=0; v2<n_gen; v2++) {
	    rf[j] += nrec(v+1,v2+1,phase[j]) * exp(gamma[v][v2] - s);
	  }
	}
      }

    } /* loop over individuals */

    /* rescale */
    for(j=0; j<n_mar-1; j++) {
      rf[j] /= (double)n_ind;
      if(rf[j] < tol/100.0) rf[j] = tol/100.0;
      else if(rf[j] > 0.5-tol/100.0) rf[j] = 0.5-tol/100.0;
    }
    /*
    Rprintf("\n"); 
    Rprintf("it: %4d ", it+1);
    for(j=0; j<n_mar-1; j++) Rprintf("%.7lf ", rf[j], "\n");
    */  

    /* check convergence */
    for(j=0, flag=0; j<n_mar-1; j++) {
      if(fabs(rf[j] - cur_rf[j]) > tol*(cur_rf[j]+tol*100.0)) {
	flag = 1; 
	break;
      }
    }

    if(!flag) break;

  } /* end EM algorithm */
  
  /*if(flag) warning("Didn't converge!\n");*/

  /* calculate log likelihood */
  *loglik = 0.0;
  for(i=0; i<n_ind; i++) { /* i = individual */
    /* initialize alpha */
    for(v=0; v<n_gen; v++) 
      alpha[v][0] = initf(v+1) + emitf(Geno[0][i], v+1, error_prob, type[0]);
    /* forward equations */
    for(j=1; j<n_mar; j++) {
      for(v=0; v<n_gen; v++) {
	alpha[v][j] = alpha[0][j-1] + 
	  stepf(1, v+1, phase[j-1], rf[j-1]);
	for(v2=1; v2<n_gen; v2++) 
	  alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
			       stepf(v2+1,v+1,phase[j-1],rf[j-1]));
	alpha[v][j] += emitf(Geno[j][i],v+1,error_prob,type[j]);
      }

    }

    curloglik = alpha[0][n_mar-1];
    for(v=1; v<n_gen; v++) 
      curloglik = addlog(curloglik, alpha[v][n_mar-1]);
    *loglik += curloglik;
  }

  if(verbose) {
    /* print final estimates */
    Rprintf(" Number of iterations to converge: %4d \n", it+1);
    for(j=0; j<n_mar-1; j++) Rprintf(" %.3lf", rf[j]);
    Rprintf("\n");
    
    Rprintf(" loglike: %10.4lf\n\n", *loglik);
  }

}












/**********************************************************************
 * 
 * calc_genoprob
 *
 * This function uses the hidden Markov model technology to calculate 
 * the genotype probabilities at each of marker and (optionally) at 
 * points between markers, conditional on all marker data for a 
 * chromosome.  This assumes data on a single chromosome
 *
 * n_ind        Number of individuals
 *
 * n_pos        Number of markers (or really positions at which to 
 *              calculate the genotype probabilities)
 *
 * n_gen        Number of different genotypes
 *  
 * geno         Genotype data, as a single vector storing the matrix 
 *              by columns, with each column corresponding to a marker
 *
 * rf           Recombination fractions
 *
 *
 * error_prob   Genotyping error probability
 *
 * genoprob     Genotype probabilities (the output); a single vector
 *              stored by columns (ind moves fastest, then mar, then
 *              genotype
 *
 * initf        Function returning log Pr(g_i)
 *
 * emitf        Function returning log Pr(O_i | g_i)
 * 
 * stepf        Function returning log Pr(g_2 | g_1)
 *
 **********************************************************************/

/* Note: true genotypes coded as 1, 2, ...
   but in the alpha's and beta's, we use 0, 1, ... */

void calc_genoprob(int n_ind, int n_pos, int *type, int *phase, int n_gen, int *geno, 
		   double *rf,  
		   double error_prob, double *genoprob, 
		   double initf(int), 
		   double emitf(int, int, double, int),
		   double stepf(int, int, int, double)) 
{
  int i, j, j2, v, v2;
  double s, **alpha, **beta;
  int **Geno;
  double ***Genoprob;
  
  /* allocate space for alpha and beta and 
     reorganize geno and genoprob */
  reorg_geno(n_ind, n_pos, geno, &Geno);
  reorg_genoprob(n_ind, n_pos, n_gen, genoprob, &Genoprob);
  allocate_alpha(n_pos, n_gen, &alpha);
  allocate_alpha(n_pos, n_gen, &beta);

  for(i=0; i<n_ind; i++) { /* i = individual */

    /* initialize alpha and beta */
    for(v=0; v<n_gen; v++) {
      alpha[v][0] = initf(v+1) + emitf(Geno[0][i], v+1, error_prob, type[0]);
      beta[v][n_pos-1] = 0.0;
    }

    /* forward-backward equations */
    for(j=1,j2=n_pos-2; j<n_pos; j++, j2--) {
      
      for(v=0; v<n_gen; v++) {
	alpha[v][j] = alpha[0][j-1] + stepf(1, v+1, phase[j-1], rf[j-1]);
	
	beta[v][j2] = beta[0][j2+1] + stepf(v+1,1, phase[j2], rf[j2]) + 
	  emitf(Geno[j2+1][i],1,error_prob, type[j2+1]);

	for(v2=1; v2<n_gen; v2++) {
	  alpha[v][j] = addlog(alpha[v][j], alpha[v2][j-1] + 
			       stepf(v2+1,v+1,phase[j-1],rf[j-1]));
	  beta[v][j2] = addlog(beta[v][j2], beta[v2][j2+1] + 
			       stepf(v+1,v2+1,phase[j2],rf[j2]) +
			       emitf(Geno[j2+1][i],v2+1,error_prob,type[j2+1]));
	}

	alpha[v][j] += emitf(Geno[j][i],v+1,error_prob,type[j]);
      }
    }

    /* calculate genotype probabilities */
    for(j=0; j<n_pos; j++) {
      s = Genoprob[0][j][i] = alpha[0][j] + beta[0][j];
      for(v=1; v<n_gen; v++) {
	Genoprob[v][j][i] = alpha[v][j] + beta[v][j];
	s = addlog(s, Genoprob[v][j][i]);
      }
      for(v=0; v<n_gen; v++) 
	Genoprob[v][j][i] = exp(Genoprob[v][j][i] - s);
    } 
  } /* loop over individuals */
  

}


void est_map_outbred(int *n_ind, int *n_mar, int *type, int *phase, int *geno, double *rf,
		     double *error_prob, double *loglik, int *maxit, 
		     double *tol, int *verbose)
{
  est_map(*n_ind, *n_mar, type, phase, 4, geno, rf, *error_prob, 
	  init_outbred, emit_outbred, step_outbred, nrec_outbred,
	  loglik, *maxit, *tol, *verbose);
}


void calc_genoprob_outbred(int *n_ind, int *n_mar, int *type,  int *phase, int *geno, 
			   double *rf, double *error_prob, double *genoprob) 

{
  calc_genoprob(*n_ind, *n_mar, type, phase, 4, geno, rf, *error_prob, genoprob,
		init_outbred, emit_outbred, step_outbred);
}     


