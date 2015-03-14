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
  File: twopts.c
  Description: Set of OneMap functions to be used with software R
               Implements the methodology of Wu et al. (2002): "Simultaneous
               maximum likelihood estimation of linkage and linkage phases in
               outcrossing species"
               Theoretical Population Biology 61, 349-363
  Written by Gabriel Rodrigues Alves Margarido
  Escola Superior de Agricultura "Luiz de Queiroz"
  Departamento de Genética - São Paulo, Brazil
  Contact: gramarga@gmail.com
  First version: 02/13/2007
  Last update: 03/02/2009
*/

/* This code was meant to be easily understood and modified, thus the matrix
     definitions were made in a very explicit way, consuming a lot of space
     in this file
   The meaning of variables and the steps of the algorithm will not be
     commented - for details check the paper cited above */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>

#define MAXITER 1000
#define THRESH log(1E-4)

/* mdrct2pt is a function used to calculate the ELEMENTWISE product of
     two matrices with dimensions 4x4 */

void mdrct2pt(double A[16], double B[16], double res[16]) {
  int i;
  for(i=0;i<16;i++) res[i] = A[i]*B[i];
}

/* mprod is a function used to calculate the PRODUCT of two matrices
     and/or vectors */

void mprod(double *A, int rowA, int colA, double *B, int rowB, int colB, double *res) {
  char *transA = "N", *transB = "N";
  int i;
  double one = 1.0, zero = 0.0;
  if(rowA > 0 && colA > 0 && rowB > 0 && colB >0) {
    F77_CALL(dgemm)(transA,transB,&rowA,&colB,&colA,&one,A,&rowA,B,&rowB,&zero,res,&rowA);
  } else
    for(i=0;i < rowA*colB;i++) res[i]=0;
}

/* THE 2 FUNCTIONS ABOVE WERE TESTED ONLY FOR THIS METHODOLOGY
   NOT TO BE USED WITH LARGE MATRICES AND/OR VECTORS */


/* H1, H2, H3 and H4 are functions that calculate the transition probability
     matrix between two markers under assignments A1, A2, A3 and A4,
     respectively, according to Wu et al. (2002) */

void H1(double r, double H[16]) {
  H[0]=(1-r)*(1-r); H[4]=r*(1-r); H[8]=H[4];  H[12]=r*r;
  H[1]=H[4];        H[5]=H[0];    H[9]=H[12]; H[13]=H[4];
  H[2]=H[4];        H[6]=H[12];   H[10]=H[0]; H[14]=H[4];
  H[3]=H[12];       H[7]=H[4];    H[11]=H[4]; H[15]=H[0];
  
  /*H[0]=(1-r)*(1-r); H[4]=r*(1-r); H[8]=r*(1-r); H[12]=r*r;
    H[1]=r*(1-r); H[5]=(1-r)*(1-r); H[9]=r*r; H[13]=r*(1-r);
    H[2]=r*(1-r); H[6]=r*r; H[10]=(1-r)*(1-r); H[14]=r*(1-r);
    H[3]=r*r; H[7]=r*(1-r); H[11]=r*(1-r); H[15]=(1-r)*(1-r);*/
}

void H2(double r, double H[16]) {
  H[0]=r*(1-r); H[4]=(1-r)*(1-r); H[8]=r*r;   H[12]=H[0];
  H[1]=H[4];    H[5]=H[0];        H[9]=H[0];  H[13]=H[8];
  H[2]=H[8];    H[6]=H[0];        H[10]=H[0]; H[14]=H[4];
  H[3]=H[0];    H[7]=H[8];        H[11]=H[4]; H[15]=H[0];
  
  /*H[0]=r*(1-r); H[4]=(1-r)*(1-r); H[8]=r*r; H[12]=r*(1-r);
    H[1]=(1-r)*(1-r); H[5]=r*(1-r); H[9]=r*(1-r); H[13]=r*r;
    H[2]=r*r; H[6]=r*(1-r); H[10]=r*(1-r); H[14]=(1-r)*(1-r);
    H[3]=r*(1-r); H[7]=r*r; H[11]=(1-r)*(1-r); H[15]=r*(1-r);*/
}

void H3(double r, double H[16]) {
  H[0]=r*(1-r); H[4]=r*r;  H[8]=(1-r)*(1-r); H[12]=H[0];
  H[1]=H[4];    H[5]=H[0]; H[9]=H[0];        H[13]=H[8];
  H[2]=H[8];    H[6]=H[0]; H[10]=H[0];       H[14]=H[4];
  H[3]=H[0];    H[7]=H[8]; H[11]=H[4];       H[15]=H[0];
  
  /*H[0]=r*(1-r); H[4]=r*r; H[8]=(1-r)*(1-r); H[12]=r*(1-r);
    H[1]=r*r; H[5]=r*(1-r); H[9]=r*(1-r); H[13]=(1-r)*(1-r);
    H[2]=(1-r)*(1-r); H[6]=r*(1-r); H[10]=r*(1-r); H[14]=r*r;
    H[3]=r*(1-r); H[7]=(1-r)*(1-r); H[11]=r*r; H[15]=r*(1-r);*/
}

void H4(double r, double H[16]) {
  H[0]=r*r;   H[4]=r*(1-r); H[8]=H[4];  H[12]=(1-r)*(1-r);
  H[1]=H[4];  H[5]=H[0];    H[9]=H[12]; H[13]=H[4];
  H[2]=H[4];  H[6]=H[12];   H[10]=H[0]; H[14]=H[4];
  H[3]=H[12]; H[7]=H[4];    H[11]=H[4]; H[15]=H[0];
  
  /*H[0]=r*r; H[4]=r*(1-r); H[8]=r*(1-r); H[12]=(1-r)*(1-r);
    H[1]=r*(1-r); H[5]=r*r; H[9]=(1-r)*(1-r); H[13]=r*(1-r);
    H[2]=r*(1-r); H[6]=(1-r)*(1-r); H[10]=r*r; H[14]=r*(1-r);
    H[3]=(1-r)*(1-r); H[7]=r*(1-r); H[11]=r*(1-r); H[15]=r*r;*/
}


/* absol is a simple function to return the absolute value of a variable */

/* double absol(double x) {
    if (x<0) return -x;
    else return x;
} */

/* if I want (a+b) and only have x=log(a) and y=log(b), log_add calculates log(a+b), that is, log[exp(x) + exp(y)] */

double log_add(double x, double y) {
  if(y > x + 300) return(y);
  else if (x > y + 300) return(x);
  if (x<y) { double t=x; x=y; y=t; }
  return (x + log1p(exp(y-x)));
}

/* if I want (a-b) and only have x=log(a) and y=log(b), log_sub calculates log(a-b), that is, log[exp(x) - exp(y)] */

double log_sub(double x, double y) {
  if(y > x + 300) return(y);
  else if (x > y + 300) return(x);
  if (x<y) { double t=x; x=y; y=t; }
  return (x + log1p(-exp(y-x)));
}

/* rf_2pt calculates the recombination fraction and log-likelihood for a given assignment */

void rf_2pt(double *I1, int p1, double *I2, int p2, int *n, int ntot, void (*Hcall)(double, double [16]), double D[16], double *rf_assign, double *log_like_assign) {
  double diff, rf, sum, *temp, *P, *num, H[16], mid[16], log_Lold, log_Lnew, two_n;
  int i, iter;
  
  num = (double *)R_alloc(p1*p2, sizeof(double));
  temp = (double *)R_alloc(p1*4, sizeof(double));
  P = (double *)R_alloc(p1*p2, sizeof(double));

  two_n = 2*ntot;

  rf = 0.25; /* the initial value for the recombination fraction */
  Hcall(rf,H);
  mprod(I1,p1,4,(double *)H,4,4,temp);
  mprod(temp,p1,4,I2,4,p2,P);
  /* the likelihood for the initial values */
  log_Lnew = 0.0;
  for (i=0;i<(p1*p2);i++)
    log_Lnew += n[i]*log(P[i]);
  diff = 1;
  
  iter = 1;
  
  /* the recombination fraction will be calculated until convergence */
  while(rf != 0.0 && rf != 1.0 && diff>THRESH && iter<MAXITER) {
    iter++;
    log_Lold = log_Lnew;
    
    mdrct2pt(D,H,mid);
    mprod(I1,p1,4,(double *)mid,4,4,temp);
    mprod(temp,p1,4,I2,4,p2,num);
    sum = 0.0;
    for (i=0;i<(p1*p2);i++) sum += n[i]*(num[i]/P[i]);
    rf = sum/(two_n); /* the new value for the recombination fraction */
    
    Hcall(rf,H);
    mprod(I1,p1,4,(double *)H,4,4,temp);
    mprod(temp,p1,4,I2,4,p2,P);
    log_Lnew = 0.0;
    for (i=0;i<(p1*p2);i++)
      log_Lnew += n[i]*log(P[i]);
    diff = log_sub(log_Lnew,log_Lold) - log_Lnew;
  }
  if(rf == 0.0) rf = 10E-10;
  if(rf == 1.0) rf = 1.0 - 10E-10;
  Hcall(rf,H);
  mprod(I1,p1,4,(double *)H,4,4,temp);
  mprod(temp,p1,4,I2,4,p2,P);
  log_Lnew = 0.0;
  for (i=0;i<(p1*p2);i++)
    log_Lnew += n[i]*log(P[i]);
	    
  *rf_assign = rf;
  *log_like_assign = log_Lnew;
}


/* r2pts is the main function of the two-point analysis proposed by this
     methodology
   it handles all the calculations under all 4 assignments */

void r2pts(double *I1, int *p1, double *I2, int *p2, int *n, int *ntot, double *r, double *log_like, double *posterior, double *LOD) {
  double D[16], H[16], *temp, *P, log_sum_like, log_like_null;
  int i;
  
  temp = (double *)R_alloc(*p1 * 4, sizeof(double));
  P = (double *)R_alloc(*p1 * *p2, sizeof(double));
  
  /* assignment 1 */
  /* this matrix has been taken directly from the paper */
  D[0]=0; D[4]=1; D[8]=1; D[12]=2;
  D[1]=1; D[5]=0; D[9]=2; D[13]=1;
  D[2]=1; D[6]=2; D[10]=0; D[14]=1;
  D[3]=2; D[7]=1; D[11]=1; D[15]=0;
  
  rf_2pt(I1,*p1,I2,*p2,n,*ntot,H1,D,&r[0],&log_like[0]);
  
  /* assignment 2 */
  D[0]=1; D[4]=0; D[8]=2; D[12]=1;
  D[1]=0; D[5]=1; D[9]=1; D[13]=2;
  D[2]=2; D[6]=1; D[10]=1; D[14]=0;
  D[3]=1; D[7]=2; D[11]=0; D[15]=1;
  
  rf_2pt(I1,*p1,I2,*p2,n,*ntot,H2,D,&r[1],&log_like[1]);
  
  /* assignment 3 */
  D[0]=1; D[4]=2; D[8]=0; D[12]=1;
  D[1]=2; D[5]=1; D[9]=1; D[13]=0;
  D[2]=0; D[6]=1; D[10]=1; D[14]=2;
  D[3]=1; D[7]=0; D[11]=2; D[15]=1;
  
  rf_2pt(I1,*p1,I2,*p2,n,*ntot,H3,D,&r[2],&log_like[2]);
  
  /* assignment 4 */
  D[0]=2; D[4]=1; D[8]=1; D[12]=0;
  D[1]=1; D[5]=2; D[9]=0; D[13]=1;
  D[2]=1; D[6]=0; D[10]=2; D[14]=1;
  D[3]=0; D[7]=1; D[11]=1; D[15]=2;
  
  rf_2pt(I1,*p1,I2,*p2,n,*ntot,H4,D,&r[3],&log_like[3]);
  
  /* the following statements calculate the posterior probability of each assignment */
  log_sum_like = -1e50;
  for (i=0;i<4;i++) log_sum_like = log_add(log_sum_like,log_like[i]);
  for (i=0;i<4;i++) posterior[i] = exp(log_like[i]-log_sum_like);
  
  /* the likelihood under null hipothesis of no linkage is calculated below */
  H1(0.5,H);
  mprod(I1,*p1,4,(double *)H,4,4,temp);
  mprod(temp,*p1,4,I2,4,*p2,P);
  log_like_null = 0.0;
  for (i=0;i<(*p1 * *p2);i++)
    log_like_null += n[i]*log(P[i]);
    
  /* the LOD score is calculated for each assignment */
  for (i=0;i<4;i++) LOD[i] = (log_like[i]-log_like_null)/log(10);
}


