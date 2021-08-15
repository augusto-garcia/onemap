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
 File: f2_est.cpp
 
 Description: Set of functions to compute the recombination fraction
 in F2 experimental populations. These functions contain the EM
 algorithms for all possible combination of types of markers
 (co-dominant and dominant).
 
 Written by Marcelo Mollinari
 
 Escola Superior de Agricultura "Luiz de Queiroz"
 Departamento de Genética - São Paulo, Brazil
 Contact: mmollina@usp.br
 First version: 09/2015
 Last update: 10/2015
 */

#include <Rcpp.h>
#include <R_ext/PrtUtil.h>
using namespace Rcpp;
using namespace std;
#define TOL 1e-6
#define LN_75 -0.28768207245178
#define rf_TOL_min 1e-50
#define rf_TOL_max 0.5-1e-50

Rcpp::NumericVector est_rf_C_C(std::vector<int> k_sub,
                               std::vector<int> k1_sub,
                               int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n1=0, n2=0, n3=0, n4=0;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
  {
    if(k_sub[k]==1)
    {
      if (k1_sub[k]==1) n1++;
      else if (k1_sub[k]==2) n2++;
      else if (k1_sub[k]==3) n3++;
      else n0++;
      
    }
    else if(k_sub[k]==2)
    {
      if (k1_sub[k]==1) n2++;
      else if (k1_sub[k]==2) n4++;
      else if (k1_sub[k]==3) n2++;
      else n0++;
    }
    else if(k_sub[k]==3)
    {
      if (k1_sub[k]==1) n3++;
      else if (k1_sub[k]==2) n2++;
      else if (k1_sub[k]==3) n1++;
      else n0++;
    }
    else n0++;
  }
  //EM algorithm
  while(abs(rold-rnew) > TOL)
  {
    rold=rnew;
    rnew=(n2+2*(n3+n4*rold*rold/((1-rold)*(1-rold)+rold*rold)))/(2*(n_ind-n0));
  }
  //Likelihood
  if(rnew > rf_TOL_max)
  {
    r(0)=rf_TOL_max;
    r(1)=0.0;
    return(r);
  }
  if(rnew < rf_TOL_min) rnew=rf_TOL_min;
  l=n4*log(rnew*rnew+(1-rnew)*(1-rnew))+
    2*n3*log(rnew)+n2*log(rnew*(1-rnew))+2*n1*log(1-rnew);
  //Likelihood unde H0: r=0.5
  l0=-M_LN2*(n4+2.0*n3+2.0*n2+2*n1);
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return(r);
}

Rcpp::NumericVector est_rf_C_D_43(std::vector<int> k_sub,
                                  std::vector<int> k1_sub,
                                  int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n1=0, n3=0, n4=0, n5=0, n6=0, n8=0;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
  {
    if(k_sub[k]==1)
    {
      if (k1_sub[k]==3) n3++;
      else if (k1_sub[k]==4) n4++;
      else n0++;
    }
    else if(k_sub[k]==2)
    {
      if (k1_sub[k]==3) n6++;
      else if (k1_sub[k]==4) n8++;
      else n0++;
    }
    else if(k_sub[k]==3)
    {
      if (k1_sub[k]==3) n1++;
      else if (k1_sub[k]==4) n5++;
      else n0++;
    }
    else n0++;
  }
  //EM algorithm
  while(abs(rold-rnew) > TOL)
  {
    rold=rnew;
    rnew=((n5 + n1 + 2*n8 + 2*n6 + n4 + n3)*rold + n1 + n4)/(2*(n_ind-n0));
  }
  //Likelihood
  if(rnew > rf_TOL_max)
  {
    r(0)=rf_TOL_max;
    r(1)=0.0;
    return(r);
  }
  if(rnew < rf_TOL_min) rnew=rf_TOL_min;
  
  l=(n1+n4)*log(rnew) + (n5 + n3)*log(1-rnew)-log(2)*n8 - log(2)*n6;
  //Likelihood unde H0: r=0.5
  l0= -M_LN2*(n_ind-n0);
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return (r);
}

Rcpp::NumericVector est_rf_C_D_51(std::vector<int> k_sub,
                                  std::vector<int> k1_sub,
                                  int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n1=0, n3=0, n4=0, n5=0, n6=0, n8=0;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
  {
    if(k_sub[k]==1)
    {
      if (k1_sub[k]==1) n1++;
      else if (k1_sub[k]==5) n5++;
      else n0++;
    }
    else if(k_sub[k]==2)
    {
      if (k1_sub[k]==1) n6++;
      else if (k1_sub[k]==5) n8++;
      else n0++;
    }
    else if(k_sub[k]==3)
    {
      if (k1_sub[k]==1) n3++;
      else if (k1_sub[k]==5) n4++;
      else n0++;
    }
    
    else n0++;
  }
  //EM algorithm
  while(abs(rold-rnew) > TOL)
  {
    rold=rnew;
    rnew=((n3 + n4 + 2*n6 + 2*n8 + n1 + n5)*rold + n4 + n1)/(2*(n_ind-n0));
  }
  //Likelihood
  if(rnew > rf_TOL_max)
  {
    r(0)=rf_TOL_max;
    r(1)=0.0;
    return(r);
  }
  if(rnew < rf_TOL_min) rnew=rf_TOL_min;
  l=(n4 + n1)*log(rnew) + (n3 + n5)*log(1 - rnew) - log(2)*n6 - log(2)*n8;
  //Likelihood unde H0: r=0.5
  l0= -M_LN2*(n_ind-n0);
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return (r);
}

Rcpp::NumericVector est_rf_D_D_43(std::vector<int> k_sub,
                                  std::vector<int> k1_sub,
                                  int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n1=0,n5=0, n11=0, n12=0;
  double r0, r1, r2;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
  {
    if(k_sub[k]==3)
    {
      if (k1_sub[k]==3) n1++;
      else if (k1_sub[k]==4) n5++;
      else n0++;
    }
    else if(k_sub[k]==4)
    {
      if (k1_sub[k]==3) n11++;
      else if (k1_sub[k]==4) n12++;
      else n0++;
    }
    else n0++;
  }
  //EM algorithm
  while(abs(rold-rnew) > TOL)
  {
    rold=rnew;
    r0=(1-rold)*(1-rold);
    r1=(1-rold)*rold;
    r2=rold*rold;
    rnew=((n5 + n11)*2.0*r1/(r2+2*r1) +
      n12*4*r1/(3*r0 + 4*r1 + 2*r2) +
      2*((n5 + n11)*r2/(r2+2*r1) +
      n12*2*r2/(3*r0 + 4*r1 + 2*r2)))/(2.0*(n_ind-n0));
  }
  r0=(1.0-rnew)*(1.0-rnew);
  r1=rnew*(1.0-rnew);
  r2=rnew*rnew;
  //Likelihood  
  if(rnew > rf_TOL_max)
  {
    r(0)=rf_TOL_max;
    r(1)=0.0;
    return(r);
  }
  if(rnew < rf_TOL_min) rnew=rf_TOL_min;
  l=n1 * (2.0*log(1.0-rnew)) +
    n5  *  (log(1.0-r0)) +
    n11  *  (log((1.0-r0) / 3.0)) +
    n12 * (log((r0 + 2.0) / 3.0));
  
  //Likelihood unde H0: r=0.5
  l0= -2 * M_LN2 *(n11 + n1) + LN_75*(n12 + n5);
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return (r);
}

Rcpp::NumericVector est_rf_D_D_51(std::vector<int> k_sub,
                                  std::vector<int> k1_sub,
                                  int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n1=0, n5=0, n11=0, n12=0;
  double r0, r1, r2;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
  {
    if(k_sub[k]==1)
    {
      if (k1_sub[k]==1) n1++;
      else if (k1_sub[k]==5) n5++;
      else n0++;
    }
    else if(k_sub[k]==5)
    {
      if (k1_sub[k]==1) n11++;
      else if (k1_sub[k]==5) n12++;
      else n0++;
    }
    else n0++;
  }
  //EM algorithm
  while(abs(rold-rnew) > TOL)
  {
    rold=rnew;
    r0=(1-rold)*(1-rold);
    r1=(1-rold)*rold;
    r2=rold*rold;
    rnew=((n5 + n11)*2.0*r1/(r2+2*r1) +
      n12*4*r1/(3*r0 + 4*r1 + 2*r2) +
      2*((n5 + n11)*r2/(r2+2*r1) +
      n12*2*r2/(3*r0 + 4*r1 + 2*r2)))/(2.0*(n_ind-n0));
  }
  
  //Likelihood
  if(rnew > rf_TOL_max)
  {
    r(0)=rf_TOL_max;
    r(1)=0.0;
    return(r);
  }
  if(rnew < rf_TOL_min) rnew=rf_TOL_min;
  l=n1 * (2.0*log(1.0-rnew)) +
    n5  *  (log(1.0-(1.0-rnew)*(1.0-rnew))) +
    n11  *  (log((1.0-(1.0-rnew)*(1.0-rnew)) / 3.0)) +
    n12 * (log(((1.0-rnew)*(1.0-rnew) + 2.0) / 3.0));
  //Likelihood unde H0: r=0.5
  l0= -2 * M_LN2 *(n11 +  n1) + LN_75*(n12 + n5);
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return (r);
}

Rcpp::NumericVector est_rf_D_D_43_51(std::vector<int> k_sub,
                                     std::vector<int> k1_sub,
                                     int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n3=0, n4=0, n9=0, n13=0;
  double r0, r1, r2;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
  {
    if(k_sub[k]==3)
    {
      if (k1_sub[k]==1) n3++;
      else if (k1_sub[k]==5) n4++;
      else n0++;
    }
    else if(k_sub[k]==4)
    {
      if (k1_sub[k]==1) n9++;
      else if (k1_sub[k]==5) n13++;
      else n0++;
    }
    else n0++;
  }
  //EM algorithm
  while(abs(rold-rnew) > TOL)
  {
    rold=rnew;
    r0=(1-rold)*(1-rold);
    r1=(1-rold)*rold;
    r2=rold*rold;
    rnew=((n4 + n9)*2*r1/(r0+2*r1) +
      n13*4*r1/(2*r0 + 4*r1 + 3*r2) +
      2*(n3 + n13*3*r2/(2*r0 + 4*r1 + 3*r2)))/(2.0*(n_ind-n0));
  }
  r0=(1.0-rnew)*(1.0-rnew);
  r1=rnew*(1.0-rnew);
  r2=rnew*rnew;
  //Likelihood
  if(rnew > rf_TOL_max)
  {
    r(0)=rf_TOL_max;
    r(1)=0.0;
    return(r);
  }
  if(rnew < rf_TOL_min) rnew=rf_TOL_min;
  l=n4  *  (log(1.0-r2)) +
    n9  *  (log((1.0-r2) / 3.0)) +
    n13 * (log((r2+2.0) / 3.0)) +
    n3 * (2.0*log(rnew));
  //Likelihood unde H0: r=0.5
  l0= -2 * M_LN2 *(n3 + n9) + LN_75*(n13 + n4);
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return (r);
}

Rcpp::NumericVector est_rf_A_A(std::vector<int> k_sub,
                               std::vector<int> k1_sub,
                               int n_ind)
{
  Rcpp::NumericVector r(2);
  int n0=0, n1=0, n2=0, n3=0, n4=0, n5=0, n6=0, n7=0, n8=0, n9=0, n10=0, n11=0, n12=0, n13=0;
  double r0, r1, r2;
  double rold=0, rnew=0.01, l, l0;
  for(int k=0; k < n_ind; k++)
  {
    if(k_sub[k]==1)
    {
      if (k1_sub[k]==1) n1++;
      else if (k1_sub[k]==2) n2++;
      else if (k1_sub[k]==3) n3++;
      else if (k1_sub[k]==4) n4++;
      else if (k1_sub[k]==5) n5++;
      else n0++;
    }
    else if(k_sub[k]==2)
    {
      if (k1_sub[k]==1) n6++;
      else if (k1_sub[k]==2) n7++;
      else if (k1_sub[k]==3) n6++;
      else if (k1_sub[k]==4) n8++;
      else if (k1_sub[k]==5) n8++;
      else n0++;
    }
    else if(k_sub[k]==3)
    {
      if (k1_sub[k]==1) n3++;
      else if (k1_sub[k]==2) n2++;
      else if (k1_sub[k]==3) n1++;
      else if (k1_sub[k]==4) n5++;
      else if (k1_sub[k]==5) n4++;
      else n0++;
    }
    else if(k_sub[k]==4)
    {
      if (k1_sub[k]==1) n9++;
      else if (k1_sub[k]==2) n10++;
      else if (k1_sub[k]==3) n11++;
      else if (k1_sub[k]==4) n12++;
      else if (k1_sub[k]==5) n13++;
      else n0++;
    }
    else if(k_sub[k]==5)
    {
      if (k1_sub[k]==1) n11++;
      else if (k1_sub[k]==2) n10++;
      else if (k1_sub[k]==3) n9++;
      else if (k1_sub[k]==4) n13++;
      else if (k1_sub[k]==5) n12++;
      else n0++;
    }
    else n0++;
  }
  //EM algorithm
  while(abs(rold-rnew) > TOL)
  {
    rold=rnew;
    r0=(1-rold)*(1-rold);
    r1=(1-rold)*rold;
    r2=rold*rold;
    rnew=((n2 + n6)+
      (n5 + n11)*2.0*r1/(r2+2*r1) +
      (n8 + n10)*r1/(r1+r0+r2) +
      (n4 + n9)*2*r1/(r0+2*r1) +
      n12*4*r1/(3*r0 + 4*r1 + 2*r2) +
      n13*4*r1/(2*r0 + 4*r1 + 3*r2) +
      2*(n3 +
      (n5 + n11)*r2/(r2+2*r1) +
      n7*r2/(r0+r2) +
      (n8 + n10)*r2/(r1+r0+r2) +
      n12*2*r2/(3*r0 + 4*r1 + 2*r2) +
      n13*3*r2/(2*r0 + 4*r1 + 3*r2)
      ))/(2.0*(n_ind-n0));
  }
  r0=(1.0-rnew)*(1.0-rnew);
  r1=rnew*(1.0-rnew);
  r2=rnew*rnew;
  //Likelihood
  if(rnew > rf_TOL_max)
  {
    r(0)=rf_TOL_max;
    r(1)=0.0;
    return(r);
  }
  if(rnew < rf_TOL_min) rnew=rf_TOL_min;
  l=n1 * (2.0*log(1.0-rnew)) +
    n2  *  (M_LN2 + log(rnew) + log(1.0-rnew)) +
    n6  *  (log(r1)) +
    n5  *  (log(1.0-r0)) +
    n11  *  (log((1.0-r0) / 3.0)) +
    n8  *  (log(1.0-r1)) +
    n10  *  (log((1.0-r1) / 3.0) + M_LN2) +
    n4  *  (log(1.0-r2)) +
    n9  *  (log((1.0-r2) / 3.0)) +
    n12 * (log((r0 + 2.0) / 3.0)) +
    n13 * (log((r2+2.0) / 3.0)) +
    n7 * (log(r2+r0)) +
    n3 * (2.0*log(rnew));
  //Likelihood unde H0: r=0.5
  l0= -2 * M_LN2 *(n3 + n9 + n11 + n6 + n1) - M_LN2 *(n7 + n10 + n2 ) +
    LN_75*(n13 + n12 + n4 + n8 + n5);
  r(0)=rnew;
  r(1)=(l-l0)/log(10.0); //transforming to base 10 logarithm
  return (r);
}
