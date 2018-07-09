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
  File: out_est.cpp 
  
  Description: Set of functions to compute the recombination fraction
  in outcross experimental populations. These functions contain the EM
  algorithms for all possible combination of types of markers
  (A, B1, B2, B3, C, D1 and D2). For more detail refer to Wu 2002. 

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
#define TOL 1e-06
#define LN3 1.098612288668109
#define LN4 1.38629436111989
#define LN_75 -0.28768207245178
#define rf_TOL_min 1e-5  
#define rf_TOL_max 1-1e-5 


Rcpp::NumericVector rf_A_A(Rcpp::NumericMatrix n,
			  int n_ind,
			  int mis)
{
  NumericVector r(8);
  int n1, n2, n3, n4;
  double l0, l;
  l0=-2.0*M_LN2*(n_ind-mis);  
  n1=n(4,1)+n(3,2)+n(2,3)+n(1,4);   
  n2=n(4,4)+n(3,3)+n(2,2)+n(1,1);
  n3=n(4,3)+n(3,4)+n(2,1)+n(1,2);
  n4=n(4,2)+n(3,1)+n(2,4)+n(1,3);  
  r(0)=(2.0*(n1)+n3+n4)/(2.0*(n_ind-mis));
  l=(n3+n4)*log((1-r(0))*r(0))+2.0*(n1)*log(r(0))+2.0*(n2)*log(1-r(0));
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/ 
  r(1)=(2.0*(n4)+n2+n1)/(2.0*(n_ind-mis));
  l=(n2+n1)*log((1-r(1))*r(1))+2.0*(n4)*log(r(1))+2.0*(n3)*log(1-r(1));
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));
  return(r);
}
Rcpp::NumericVector rf_A_B1(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = - (M_LN2 * (n(2,1)+n(1,1) + n(4,1)+n(3,1)+2*n(2,2)+2*n(1,3)) +
	  2 * M_LN2*(n(4,2)+n(3,3)+n(2,3)+n(1,2)+n(4,3)+n(3,2)));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;

  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(4,1)+n(3,1)+n(2,1)+n(1,1)) +
	    2*(n(2,2)+n(1,3)) + 
	    n(4,2)+n(4,1)+n(3,3)+n(3,1)+n(2,3)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=(n(4,2)+n(3,3)+n(2,3)+n(1,2))*log(rnew-(rnew*rnew))+(n(4,1)+n(3,1)+2*n(2,2)+2*n(1,3))*log(rnew)+(2*n(4,3)+2*n(3,2)+n(2,1)+n(1,1))*log(1-rnew);
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(4,1)+n(3,1)+n(2,1)+n(1,1)) + 
	    2*(n(2,3)+n(1,2))+
	    n(4,3)+n(4,1)+n(3,2)+n(3,1)+n(2,2)+n(1,3))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=(n(4,3)+n(3,2)+n(2,2)+n(1,3))*log(rnew-(rnew*rnew))+(n(4,1)+n(3,1)+2*n(2,3)+2*n(1,2))*log(rnew)+(2*n(4,2)+2*n(3,3)+n(2,1)+n(1,1))*log(1-rnew);
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_A_B2(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = - (M_LN2 * (n(2,1)+n(1,1) + n(4,1)+n(3,1)+2*n(2,2)+2*n(1,3)) +
	  2 * M_LN2*(n(4,2)+n(3,3)+n(2,3)+n(1,2)+n(4,3)+n(3,2)));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;

  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(4,1)+n(3,1)+n(2,1)+n(1,1)) +
	    2*(n(3,2)+n(1,3)) + 
	    n(4,2)+n(4,1)+n(3,3)+n(2,1)+n(2,3)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=(n(4,2)+n(3,3)+n(2,3)+n(1,2))*log(rnew-(rnew*rnew))+(n(4,1)+2*n(3,2)+n(2,1)+2*n(1,3))*log(rnew)+(2*n(4,3)+2*n(2,2)+n(3,1)+n(1,1))*log(1-rnew);
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(4,1)+n(3,1)+n(2,1)+n(1,1)) + 
	    2*(n(2,3)+n(4,2))+
	    n(4,3)+n(3,2)+n(3,1)+n(2,2)+n(1,3)+n(1,1))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=(n(4,3)+n(3,2)+n(2,2)+n(1,3))*log(rnew-(rnew*rnew))+(2*n(4,2)+n(3,1)+2*n(2,3)+n(1,1))*log(rnew)+(n(4,1)+2*n(3,3)+n(2,1)+2*n(1,2))*log(1-rnew);
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_A_B3(Rcpp::NumericMatrix n,
			    int n_ind,
			    int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0= -2.0*M_LN2*(n_ind-mis);

  /*EM algorithm*/
  rold=0, rnew=0.01;

  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(((n(3,2)+n(2,2))*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 	    
	    2.0*(n(1,3)+n(4,1))+n(3,3)+n(2,3)+n(4,2)+n(1,2)+n(3,1)+n(2,1))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=(n(3,2)+n(2,2))*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
    (n(3,3)+n(2,3)+n(4,2)+n(1,2)+n(3,1)+n(2,1))*log(rnew-(rnew*rnew))+
    2.0*(n(1,3)+n(4,1))*log(rnew)+
    2.0*(n(4,3)+n(1,1))*log(1.0-rnew);
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(((n(4,2)+n(1,2))*(rold*rold))/((rold*rold)/2.0+((1-rold)*(1-rold))/2.0) + 
	    2.0*(n(2,3)+n(3,1)) + n(4,3)+n(1,3)+n(3,2)+n(2,2)+n(4,1)+n(1,1))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=(n(4,2)+n(1,2))*log((2.0*(rnew*rnew)-2.0*rnew+1)/2.0)+
    (n(4,3)+n(1,3)+n(3,2)+n(2,2)+n(4,1)+n(1,1))*log(rnew-(rnew*rnew))+
    2.0*(n(2,3)+n(3,1))*log(rnew)+
    2.0*(n(3,3)+n(2,1))*log(1.0-rnew);
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/

  r(2)=abs(1.0-r(1));

  r(3)=abs(1.0-r(0));  

  return(r);
}

Rcpp::NumericVector rf_A_C(Rcpp::NumericMatrix n,
			   int n_ind,
			   int mis)
{
  NumericVector r(8);
  double l, l0, r0, r1, r2, rnew, rold;
  //Likelihoods under h0: r=0.5
  l0=-2*M_LN2*n(4,2)-(LN4-LN3)*n(4,1)-2*M_LN2*(n(3,2)+n(2,2))-(LN4-LN3)*(n(3,1)+n(2,1))-2*M_LN2*n(1,2)-(LN4-LN3)*n(1,1);
  //EM algorithm
  rold=0, rnew=0.01;
  
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      r0=(1-rold)*(1-rold);
      r1=(1-rold)*rold;
      r2=rold*rold;
      rnew=((2.0*n(1,1)*r1)/(2.0*r1+r0) +
	    (n(4,1)*(2.0*r2+2.0*r1))/(r2+2.0*r1) +
	    ((n(3,1)+n(2,1))*(2.0*r2+r1))/(r2+r1+r0) + 
	    n(3,2)+n(2,2)+2.0*n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=(n(3,1)+n(2,1))*log(rnew*rnew-rnew+1)+n(4,1)*log(2.0*rnew-rnew*rnew)+(n(3,2)+n(2,2))*log(rnew-rnew*rnew)+n(1,1)*log(1-rnew*rnew)+2.0*n(1,2)*log(rnew)+2.0*n(4,2)*log(1-rnew);
  r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      r0=(1-rold)*(1-rold);
      r1=(1-rold)*rold;
      r2=rold*rold;
      rnew=((2.0*n(2,1)*r1)/(2.0*r1+r0) +
	    (n(3,1)*(2.0*r2+2.0*r1))/(r2+2.0*r1) +
	    ((n(4,1)+n(1,1))*(2.0*r2+r1))/(r2+r1+r0) + 
	    n(4,2)+2.0*n(2,2)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=(n(4,1)+n(1,1))*log(rnew*rnew-rnew+1)+n(3,1)*log(2.0*rnew-rnew*rnew)+(n(4,2)+n(1,2))*log(rnew-rnew*rnew)+n(2,1)*log(1-rnew*rnew)+2.0*n(2,2)*log(rnew)+2.0*n(3,2)*log(1-rnew);
  r(5)=r(6)=(l-l0)/log(10.0); //transforming to base 10 logarithm

  r(2)=abs(1.0-r(1));

  r(3)=abs(1.0-r(0));  

  return(r);
} 
Rcpp::NumericVector rf_A_D1(Rcpp::NumericMatrix n,
			   int n_ind,
			   int mis)
{
  NumericVector r(8);
  double l, l0,  rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0, rnew=0.01;
  
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(4,2)+n(3,2)+n(2,1)+n(1,1))+(rold+1)*(n(4,1)+n(3,1)+n(2,2)+n(1,2)))/(2.0*(n_ind-mis));
    }
  r(1)=r(0)=rnew;
  l=(n(4,1)+n(3,1)+n(2,2)+n(1,2))*log(rnew)+(n(4,2)+n(3,2)+n(2,1)+n(1,1))*log(1-rnew);
  r(4)=r(5)=r(6)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(3)=r(2)=abs(1.0-r(1));
  return(r);
}
Rcpp::NumericVector rf_A_D2(Rcpp::NumericMatrix n,
			   int n_ind,
			   int mis)
{
  NumericVector r(8);
  double l, l0,  rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0, rnew=0.01;
  
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(4,2)+n(3,1)+n(2,2)+n(1,1))+(rold+1)*(n(4,1)+n(3,2)+n(2,1)+n(1,2)))/(2.0*(n_ind-mis));
    }
  r(2)=r(0)=rnew;
  l=(n(4,1)+n(3,2)+n(2,1)+n(1,2))*log(rnew)+(n(4,2)+n(3,1)+n(2,2)+n(1,1))*log(1-rnew);
  r(4)=r(5)=r(6)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(3)=r(1)=abs(1.0-r(0));
  return(r);
}
Rcpp::NumericVector rf_B1_B1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(2.0*n(3,3)+2.0*n(3,2)+n(3,1)+2.0*n(2,3)+2.0*n(2,2)+n(2,1)+2.0*n(1,3)+2.0*n(1,2)+n(1,1));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,1)+n(2,1)+n(1,3)+n(1,2)+n(1,1))*rold+n(3,2)+n(3,1)+n(2,3)+n(2,1)+n(1,3)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=((n(3,2)+n(2,3))*log(rnew-rnew*rnew)+(n(3,1)+n(2,1))*log(rnew)+(n(1,3)+n(1,2))*log(rnew/2.0)+(2.0*n(3,3)+2.0*n(2,2)+n(1,1))*log(1-rnew));
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,1)+n(2,1)+n(1,3)+n(1,2)+n(1,1))*rold+n(3,3)+n(3,1)+n(2,2)+n(2,1)+n(1,3)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=((n(3,3)+n(2,2))*log(rnew-rnew*rnew)+(n(3,1)+n(2,1))*log(rnew)+(n(1,3)+n(1,2))*log(rnew/2.0)+(2.0*n(3,2)+2.0*n(2,3)+n(1,1))*log(1-rnew));
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_B1_B2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(2*n(3,3)+n(2,1)+n(3,1)+2*n(2,2)+n(1,1))-2*M_LN2*(n(3,2)+n(2,3)+n(1,3)+n(1,2));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,1)+n(2,1)+n(1,3)+n(1,2)+2*n(1,1))*rold+n(3,2)+n(3,1)+n(2,3)+2.0*n(2,2)+n(1,3))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=(n(3,2)+n(2,3))*log(rnew-rnew*rnew)+(n(3,1)+2*n(2,2))*log(rnew)+n(1,3)*log(rnew/2.0)+n(1,2)*log(-(rnew-1.0)/2.0)+(2.0*n(3,3)+n(2,1))*log(1-rnew)-M_LN2*n(1,1);
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,1)+n(2,1)+n(1,3)+n(1,2)+2*n(1,1))*rold+n(3,3)+2*n(3,2)+n(2,2)+n(2,1)+n(1,3))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=(n(3,3)+n(2,2))*log(rnew-rnew*rnew)+(n(2,1)+2*n(3,2))*log(rnew)+n(1,3)*log(rnew/2.0)+n(1,2)*log(-(rnew-1.0)/2.0)+(2.0*n(2,3)+n(3,1))*log(1-rnew)-M_LN2*n(1,1);
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_B1_B3(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(2.0*n(3,3)+n(3,2)+2.0*n(3,1)+2.0*n(2,3)+n(2,2)+2.0*n(2,1)+2.0*n(1,3)+n(1,2)+2.0*n(1,1));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((2.0*n(1,3)+4*n(1,2)+2.0*n(1,1))*(rold*rold*rold) +
	    (2.0*n(3,2)+4*n(3,1)+2.0*n(2,3)+2.0*n(2,2)+2.0*n(2,1)-4*n(1,2)-2.0*n(1,1))*(rold*rold) +
	    (-2.0*n(3,2)-4*n(3,1)-2.0*n(2,3)-2.0*n(2,1)-n(1,3)+2.0*n(1,2)+n(1,1))*rold + 
	    n(3,2)+2.0*n(3,1)+n(2,3)+n(2,1)+n(1,3))/((2.0*rold*rold-2.0*rold+1)*2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=n(2,2)*log(2.0*rnew*rnew-2.0*rnew+1)+(n(2,3)+n(2,1))*log(rnew-rnew*rnew)+n(3,2)*log(2.0*rnew-2.0*rnew*rnew)+2.0*n(3,1)*log(rnew)+n(1,3)*log(rnew/2.0)+n(1,1)*log(-(rnew-1)/2.0)+2.0*n(3,3)*log(1-rnew)-M_LN2*n(1,2);
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	     


  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((2.0*n(1,3)+4*n(1,2)+2.0*n(1,1))*(rold*rold*rold) +
	    (2.0*n(3,3)+2.0*n(3,2)+2.0*n(3,1)+2.0*n(2,2)+4*n(2,1)-4*n(1,2)-2.0*n(1,1))*(rold*rold) +
	    (-2*n(3,3)-2*n(3,1)-2*n(2,2)-4*n(2,1)-n(1,3)+2*n(1,2)+n(1,1))*rold + 
	    +n(3,3)+n(3,1)+n(2,2)+2*n(2,1)+n(1,3))/((2.0*rold*rold-2.0*rold+1)*2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=n(3,2)*log(2.0*rnew*rnew-2.0*rnew+1)+(n(3,3)+n(3,1))*log(rnew-rnew*rnew)+n(2,2)*log(2.0*rnew-2.0*rnew*rnew)+2.0*n(2,1)*log(rnew)+n(1,3)*log(rnew/2.0)+n(1,1)*log(-(rnew-1)/2.0)+2.0*n(2,3)*log(1-rnew)-M_LN2*n(1,2);
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_B1_C(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = LN_75*(n(3,1) + n(2,1) + n(1,1)) - 2.0*M_LN2*(n(3,2)+n(2,2)+n(1,2));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(1,2)+n(1,1))*rold*rold*rold*rold+
	    (n(2,2)+n(2,1)-2*n(1,2)-4*n(1,1))*rold*rold*rold+
	    (-2*n(3,1)-3*n(2,2)-n(2,1)+4*n(1,1))*rold*rold+
	    (2*n(3,1)+3*n(2,2)-2*n(2,1)+n(1,2)-3*n(1,1))*rold-2*n(3,1)-2*n(2,2)-2*n(1,2))/
	(((rold-2)*(rold*rold-rold+1))*2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=n(2,1)*log(rnew*rnew-rnew+1)+
    n(3,1)*log(2*rnew-rnew*rnew)+
    n(2,2)*log(rnew-rnew*rnew)+
    n(1,2)*log(rnew/2)+
    n(1,1)*log(-(rnew-2)/2)+
    2*n(3,2)*log(1-rnew);
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	       
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(1,2)+n(1,1))*rold*rold*rold*rold+
	    (n(3,2)+n(3,1)-2*n(1,2)-4*n(1,1))*rold*rold*rold+
	    (-3*n(3,2)-n(3,1)-2*n(2,1)+4*n(1,1))*rold*rold+
	    (3*n(3,2)-2*n(3,1)+2*n(2,1)+n(1,2)-3*n(1,1))*rold-2*n(3,2)-2*n(2,1)-2*n(1,2))/
	(((rold-2)*(rold*rold-rold+1))*2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=n(3,1)*log(rnew*rnew-rnew+1)+
    n(2,1)*log(2*rnew-rnew*rnew)+
    n(3,2)*log(rnew-rnew*rnew)+
    n(1,2)*log(rnew/2)+
    n(1,1)*log(-(rnew-2)/2)+
    2*n(2,2)*log(1-rnew);
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_B1_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(3,2)+n(3,1)+n(2,2)+n(2,1)+n(1,2)+n(1,1))+
	    n(3,1)+n(2,1)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=r(1)=rnew;
  r(2)=r(3)=abs(1.0-r(0));
  l=(n(3,1)+n(2,1)+n(1,2))*log(rnew)+(n(3,2)+n(2,2)+n(1,1))*log(1-rnew);
  r(5)=r(6)=r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  return(r);
}
Rcpp::NumericVector rf_B1_D2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(rold*(n(3,2)+n(3,1)+n(2,2)+n(2,1)+2*(n(1,2)+n(1,1)))+
		 n(3,1)+n(2,2))/(2.0*(n_ind-mis));
    }
  r(0)=r(2)=rnew;
  r(1)=r(3)=abs(1.0-r(0));
  l=(n(3,1)+n(2,2))*log(rnew)+(n(3,2)+n(2,1))*log(1-rnew)-M_LN2*(n(1,2) + n(1,1));
  r(5)=r(6)=r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  return(r);
}
Rcpp::NumericVector rf_B2_B2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(2.0*n(3,3)+2.0*n(3,2)+n(3,1)+2.0*n(2,3)+2.0*n(2,2)+n(2,1)+2.0*n(1,3)+2.0*n(1,2)+n(1,1));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,1)+n(2,1)+n(1,3)+n(1,2)+n(1,1))*rold+
	    n(3,2)+n(3,1)+n(2,3)+n(2,1)+n(1,3)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=(n(3,2)+n(2,3))*log(rnew-rnew*rnew)+
    (n(3,1)+n(2,1))*log(rnew)+(n(1,3)+n(1,2))*log(rnew/2.0)+
    (2.0*n(3,3)+2.0*n(2,2)+n(1,1))*log(1-rnew);
  r(4)=r(7)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  rold=0, rnew=0.01;	     
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,1)+n(2,1)+n(1,3)+n(1,2)+n(1,1))*rold+2.0*(n(3,2)+n(2,3))+
	    n(3,3)+n(2,2)+n(1,1))/(2.0*(n_ind-mis));
    }
  r(1)=rnew;
  l=((n(3,3)+n(2,2))*log(rnew-rnew*rnew)+(2.0*(n(3,2)+n(2,3))+n(1,1))*log(rnew)+(n(1,3)+n(1,2))*log(-(rnew-1)/2)+(n(3,1)+n(2,1))*log(1-rnew));
  r(5)=r(6)=(l-l0)/log(10.0); /*transforming to base 10 logarithm*/
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_B2_B3(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(2.0*n(3,3)+n(3,2)+2.0*n(3,1)+2.0*n(2,3)+n(2,2)+2.0*n(2,1)+2.0*n(1,3)+n(1,2)+2.0*n(1,1));
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((2.0*n(1,3)+4*n(1,2)+2.0*n(1,1))*rold*rold*rold+
	    (2.0*n(3,2)+4*n(3,1)+2.0*n(2,3)+2.0*n(2,2)+2.0*n(2,1)-4*n(1,2)-2.0*n(1,1))*rold*rold+
	    (-2.0*n(3,2)-4*n(3,1)-2.0*n(2,3)-2.0*n(2,1)-n(1,3)+2.0*n(1,2)+n(1,1))*rold+
	    n(3,2)+2.0*n(3,1)+n(2,3)+n(2,1)+n(1,3))/
	(2.0*(n_ind-mis)*(2.0*(rold*rold)-2.0*rold+1));
    }
  r(0)=rnew;
  l=n(2,2)*log(2*rnew*rnew-2*rnew+1)+(n(2,3)+n(2,1))*log(rnew-rnew*rnew)+n(3,2)*log(2*rnew-2*rnew*rnew)+2*n(3,1)*log(rnew)+n(1,3)*log(rnew/2)+n(1,1)*log(-(rnew-1)/2)+2*n(3,3)*log(1-rnew)-M_LN2*n(1,2);

  r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  rold=0, rnew=0.01;	    
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((2*n(1,3)+4*n(1,2)+2*n(1,1))*rold*rold*rold+
	    (2*n(3,3)+2*n(3,2)+2*n(3,1)+4*n(2,3)+2*n(2,2)-2*n(1,3)-4*n(1,2))*rold*rold+
	    (-2*n(3,3)-2*n(3,1)-4*n(2,3)-2*n(2,2)+n(1,3)+2*n(1,2)-n(1,1))*rold+
	    n(3,3)+n(3,1)+2*n(2,3)+n(2,2)+n(1,1))/
	(2.0*(n_ind-mis)*(2.0*(rold*rold)-2.0*rold+1));
    }
  r(1)=rnew;
  l=n(3,2)*log(2*(rnew*rnew)-2*rnew+1)+
    (n(3,3)+n(3,1))*log(rnew-(rnew*rnew))+
    n(2,2)*log(2*rnew-2*(rnew*rnew))+
    2*n(2,3)*log(rnew)+n(1,1)*log(rnew/2)+
    n(1,3)*log(-(rnew-1)/2)+
    2*n(2,1)*log(1-rnew)-
    M_LN2*n(1,2);
  r(5)=r(6)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  

  return(r);
}
Rcpp::NumericVector rf_B2_C(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  //Likelihoods under h0: r=0.5
  l0 = LN_75*(n(3,1) + n(2,1) + n(1,1)) - 2.0*M_LN2*(n(3,2) + n(2,2) + n(1,2));
  //EM algorithm
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(1,2)+n(1,1))*rold*rold*rold*rold+
	    (n(2,2)+n(2,1)-2.0*n(1,2)-4*n(1,1))*rold*rold*rold+
	    (-2.0*n(3,1)-3*n(2,2)-n(2,1)+4*n(1,1))*rold*rold+
	    (2.0*n(3,1)+3*n(2,2)-2.0*n(2,1)+n(1,2)-3*n(1,1))*rold-
	    2.0*n(3,1)-2.0*n(2,2)-2.0*n(1,2))/(2.0*(n_ind-mis) * (rold-2.0)*((rold*rold)-rold+1));
    }
  r(0)=rnew;
  l=n(2,1)*log((rnew*rnew)-rnew+1)+n(3,1)*log(2.0*rnew-(rnew*rnew))+
    n(2,2)*log(rnew-(rnew*rnew))+n(1,2)*log(rnew/2.0)+n(1,1)*log(-(rnew-2)/2)+
    2.0*n(3,2)*log(1-rnew);
    r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  rold=0, rnew=0.01;	       
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(1,2)+n(1,1))*rold*rold*rold*rold+
	    (n(3,2)+n(3,1)+2.0*n(2,2)+2.0*n(2,1)+2.0*n(1,1))*rold*rold*rold+
	    (2.0*n(3,1)-2.0*n(2,1)-2.0*n(1,1))*rold*rold+
	    (n(3,1)+2.0*n(2,1)+n(1,2)+3*n(1,1))*rold+
	    n(3,2)+2.0*n(2,2))/(2.0*(n_ind-mis) * (rold+1.0)*((rold*rold)-rold+1));
    }
  r(1)=rnew;
  l=n(3,1)*log(rnew*rnew-rnew+1)+
    n(3,2)*log(rnew-rnew*rnew)+
    n(2,1)*log(1-rnew*rnew)+
    n(1,1)*log((rnew+1)/2)+
    2.0*n(2,2)*log(rnew)+
    n(1,2)*log(-(rnew-1)/2.0);
  r(5)=r(6)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_B2_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,2)+n(3,1)+n(2,2)+n(2,1)+2*n(1,2)+2*n(1,1))*rold+n(3,1)+n(2,2))/(2.0*(n_ind-mis));
    }
  r(0)=r(1)=rnew;
  r(2)=r(3)=abs(1.0-r(0));
  l=(n(3,1)+n(2,2))*log(rnew)+(n(3,2)+n(2,1))*log(1-rnew)-M_LN2*(n(1,2)+n(1,1));
  r(5)=r(6)=r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return(r);
}
Rcpp::NumericVector rf_B2_D2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,2)+n(3,1)+n(2,2)+n(2,1)+n(1,2)+n(1,1))*rold+n(3,1)+n(2,1)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=r(2)=rnew;
  r(1)=r(3)=abs(1.0-r(0));
  l=(n(3,1)+n(2,1)+n(1,2))*log(rnew)+(n(3,2)+n(2,2)+n(1,1))*log(1-rnew);
  r(5)=r(6)=r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return(r);
}
Rcpp::NumericVector rf_B3_B3(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  //Likelihoods under h0: r=0.5
  l0 = -M_LN2*(2*n(3,3)+n(3,2)+2*n(3,1)+2*n(2,3)+n(2,2)+2*n(2,1)+2*n(1,3)+n(1,2)+2*n(1,1));
    //EM algorithm
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((2*n(3,2)+4*n(3,1)+2*n(2,3)+2*n(2,2)+2*n(2,1)+4*n(1,3)+2*n(1,2))*rold*rold+
	    (-2*n(3,2)-4*n(3,1)-2*n(2,3)-2*n(2,1)-4*n(1,3)-2*n(1,2))*rold+
	    n(3,2)+2*n(3,1)+n(2,3)+n(2,1)+2*n(1,3)+n(1,2))/
	(2.0*(n_ind-mis)*(2*(rold*rold)-2*rold+1));
    }
  r(0)=rnew;
  
  if(rnew < rf_TOL_min) rnew=rf_TOL_min; /*to avoid operations with zero*/
  if(rnew > rf_TOL_max) rnew=rf_TOL_max; /*to avoid operations with Inf*/
  
  l=n(2,2)*log(2*(rnew*rnew)-2*rnew+1)+(n(2,3)+n(2,1))*log(rnew-(rnew*rnew))+
    (n(3,2)+n(1,2))*log(2*rnew-2*(rnew*rnew))+(2*n(3,1)+2*n(1,3))*log(rnew)+
    (2*n(3,3)+2*n(1,1))*log(1-rnew);
    r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
    rold=0, rnew=0.01;	       
    while(abs(rold-rnew) > TOL)
      {
	rold=rnew;
	rnew=((2*n(3,3)+2*n(3,2)+2*n(3,1)+2*n(2,3)+2*n(2,2)+2*n(2,1)+2*n(1,3)+2*n(1,2)+2*n(1,1))*rold*rold+
	      (-2*n(3,3)-2*n(3,1)-2*n(2,2)-2*n(1,3)-2*n(1,1))*rold+
	      n(3,3)+n(3,1)+n(2,2)+n(1,3)+n(1,1))/
	  (2.0*(n_ind-mis)*(2*(rold*rold)-2*rold+1));      
      }
    r(2)=r(1)=rnew;
    l=(n(3,2)+n(1,2))*log(2*(rnew*rnew)-2*rnew+1)+(n(2,3)+n(2,1))*log((2*(rnew*rnew)-2*rnew+1)/2)+
      (n(3,3)+n(3,1)+n(1,3)+n(1,1))*log(rnew-(rnew*rnew))+n(2,2)*log(2*rnew-2*(rnew*rnew));
  r(5)=r(6)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_B3_C(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  //Likelihoods under h0: r=0.5
  l0 = LN_75*(n(3,1)+n(2,1)+n(1,1)) - 2*M_LN2*(n(3,2)+n(2,2)+n(1,2));
  //EM algorithm
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(2,2)+n(2,1)+2*n(1,2)+2*n(1,1))*rold*rold*rold*rold+
	    (-2*n(3,1)-2*n(2,2)-4*n(1,2)-6*n(1,1))*rold*rold*rold+
	    (6*n(1,1)-3*n(2,1))*rold*rold+
	    (n(2,2)-2*n(2,1)+2*n(1,2)-4*n(1,1))*rold-
	    2*n(3,1)-2*n(2,2)-4*n(1,2))/
	(2.0*(n_ind-mis)*((rold-2)*(rold+1)*((rold*rold)-rold+1)));
    }
  r(0)=rnew;
  l=n(2,1)*log(rnew*rnew-rnew+1)+
    n(3,1)*log(2*rnew-rnew*rnew)+
    n(2,2)*log(rnew-(rnew*rnew))+
    n(1,1)*log(1-(rnew*rnew))+
    2*n(1,2)*log(rnew)+2*n(3,2)*log(1-rnew);
  r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  rold=0, rnew=0.01;	       
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((4*n(3,2)+4*n(3,1)+4*n(2,2)+4*n(2,1)+4*n(1,2)+4*n(1,1))*rold*rold*rold*rold*rold*rold+
	    (-12*n(3,2)-4*n(3,1)-8*n(2,2)-16*n(2,1)-12*n(1,2)-4*n(1,1))*rold*rold*rold*rold*rold+
	    (16*n(3,2)-4*n(3,1)+6*n(2,2)+26*n(2,1)+16*n(1,2)-4*n(1,1))*rold*rold*rold*rold+
	    (-12*n(3,2)+4*n(3,1)-2*n(2,2)-26*n(2,1)-12*n(1,2)+4*n(1,1))*rold*rold*rold+
	    (3*n(3,2)-n(3,1)-2*n(2,2)+14*n(2,1)+3*n(1,2)-n(1,1))*rold*rold+
	    (n(3,2)-n(3,1)-4*n(2,1)+n(1,2)-n(1,1))*rold-
	    n(3,2)-n(1,2))/(2*(n_ind-mis)*((rold*rold)-rold+1)*(2*(rold*rold)-2*rold-1)*(2*(rold*rold)-2*rold+1));
    }
  r(2)=r(1)=rnew;
  l=n(2,2)*log((2*(rnew*rnew)-2*rnew+1)/2)+
    n(2,1)*log(-(2*(rnew*rnew)-2*rnew-1)/2)+
    (n(3,1)+n(1,1))*log((rnew*rnew)-rnew+1)+
    (n(3,2)+n(1,2))*log(rnew-(rnew*rnew));
  r(5)=r(6)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_B3_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,2)+n(3,1)+2*n(2,2)+2*n(2,1)+n(1,2)+n(1,1))*rold+n(3,1)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=r(1)=rnew;
  r(2)=r(3)=abs(1.0-r(0));
  l=(n(3,1)+n(1,2))*log(rnew)+(n(3,2)+n(1,1))*log(1-rnew)-log(2.0)*n(2,2)-log(2.0)*n(2,1);
  r(4)=r(5)=r(6)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return(r);
}
Rcpp::NumericVector rf_B3_D2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(3,2)+n(3,1)+2*n(2,2)+2*n(2,1)+n(1,2)+n(1,1))*rold+n(3,1)+n(1,2))/(2.0*(n_ind-mis));
    }
  r(0)=r(2)=rnew;
  r(1)=r(3)=abs(1.0-r(0));
  l=(n(3,1)+n(1,2))*log(rnew)+(n(3,2)+n(1,1))*log(1-rnew)-log(2.0)*n(2,2)-log(2.0)*n(2,1);
  r(4)=r(5)=r(6)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return(r);
}
Rcpp::NumericVector rf_C_C(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  //Likelihoods under h0: r=0.5
  l0 = LN_75*(n(2,1)+n(1,1))-2*M_LN2*(n(2,2)+n(1,2));
  //EM algorithm
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=-((n(2,1)+n(1,2)-2*n(1,1))*rold*rold+
	     (-2*n(2,1)-2*n(1,2)+4*n(1,1))*rold+
	     3*n(2,1)+3*n(1,2))/
	((n_ind-mis)*(rold-2)*((rold*rold)-2*rold+3));
    }
  r(0)=rnew;
  l=n(1,1)*log(((rnew*rnew)-2*rnew+3)/3)+
    n(1,2)*log(-((rnew*rnew)-2*rnew)/3)+
    n(2,1)*log(2*rnew-(rnew*rnew))+
    2*n(2,2)*log(1-rnew);
  r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  rold=0, rnew=0.01;	       
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(2,2)+n(2,1)+n(1,2)+n(1,1))*rold*rold*rold*rold+
	    (-2*n(2,2)-6*n(1,1))*rold*rold*rold+
	    (-3*n(2,1)-3*n(1,2)+6*n(1,1))*rold*rold+
	    (n(2,2)-2*n(2,1)-2*n(1,2)-5*n(1,1))*rold
	    -2*n(2,2))/(2.0*(n_ind-mis)*(rold-2)*(rold+1)*((rold*rold)-rold+1));
    }
  r(2)=r(1)=rnew;
  l=n(2,1)*log((rnew*rnew)-rnew+1)+
    n(1,2)*log(((rnew*rnew)-rnew+1)/3)+
    n(1,1)*log(-((rnew*rnew)-rnew-2)/3)+
    n(2,2)*log(rnew-(rnew*rnew));
  r(5)=r(6)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  r(3)=abs(1.0-r(0));  
  return(r);
}
Rcpp::NumericVector rf_C_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(2,2)+n(2,1)+n(1,2)+n(1,1))*rold*rold*rold+
	    (-n(2,2)+n(1,2)-2*n(1,1))*rold*rold+
	    (-2*n(2,2)-3*n(2,1)-6*n(1,2)-3*n(1,1))*rold-2*n(2,1))/
	(2.0*(n_ind-mis)*(rold-2)*(rold+1));
    }
  r(0)=r(1)=rnew;
  r(2)=r(3)=abs(1.0-r(0));
  l=n(1,2)*log((rnew+1)/3)+n(2,1)*log(rnew)+n(1,1)*log(-(rnew-2)/3)+n(2,2)*log(1-rnew);
  r(5)=r(6)=r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return(r);
}
Rcpp::NumericVector rf_C_D2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=((n(2,2)+n(2,1)+n(1,2)+n(1,1))*rold*rold*rold+
	    (-n(2,2)+n(1,2)-2*n(1,1))*rold*rold+
	    (-2*n(2,2)-3*n(2,1)-6*n(1,2)-3*n(1,1))*rold-2*n(2,1))/
	(2.0*(n_ind-mis)*(rold-2)*(rold+1));
    }
  r(0)=r(2)=rnew;
  r(1)=r(3)=abs(1.0-r(0));
  l=n(1,2)*log((rnew+1)/3)+n(2,1)*log(rnew)+n(1,1)*log(-(rnew-2)/3)+n(2,2)*log(1-rnew);
  r(5)=r(6)=r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return(r);
}
Rcpp::NumericVector rf_D1_D1(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(n(2,2)*rold+n(2,1)*rold+n(1,2)*rold+n(1,1)*rold+n(2,1)+n(1,2))/
	(2.0*(n_ind-mis));
    }
  r(0)=r(1)=rnew;
  r(2)=r(3)=abs(1.0-r(0));
  l=(n(2,1)+n(1,2))*log(rnew)+(n(2,2)+n(1,1))*log(1-rnew);
  r(5)=r(6)=r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return(r);
}
Rcpp::NumericVector rf_D2_D2(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  /*Likelihoods under h0: r=0.5*/
  l0 = -M_LN2*(n_ind-mis);
  /*EM algorithm*/
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(n(2,2)*rold+n(2,1)*rold+n(1,2)*rold+n(1,1)*rold+n(2,1)+n(1,2))/
	(2.0*(n_ind-mis));
    }
  r(0)=r(2)=rnew;
  r(1)=r(3)=abs(1.0-r(0));
  l=(n(2,1)+n(1,2))*log(rnew)+(n(2,2)+n(1,1))*log(1-rnew);
  r(5)=r(6)=r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return(r);
}

/*
Rcpp::NumericVector rf_X_X(Rcpp::NumericMatrix n,
			     int n_ind,
			     int mis)
{
  NumericVector r(8);
  double l, l0, rnew, rold;
  //Likelihoods under h0: r=0.5
  l0 = 
    //EM algorithm
  rold=0;
  rnew=0.01;
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(
	    
	    
	    
	    )/(2.0*(n_ind-mis));
    }
  r(0)=rnew;
  l=
    
    
    
    
    r(4)=r(7)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  rold=0, rnew=0.01;	       
  while(abs(rold-rnew) > TOL)
    {
      rold=rnew;
      rnew=(
	    
	    
	    )/(2.0*(n_ind-mis));
      
    }
  r(1)=rnew;
  l=



    r(5)=r(6)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  r(2)=abs(1.0-r(1));
  r(3)=abs(1.0-r(0));  
  return(r);
}
*/
