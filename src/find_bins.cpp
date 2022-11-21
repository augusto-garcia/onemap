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
 File: find_bins.c
 Description: Search for markers with same inforation and attributes them to bins
              The output is a vector containing the bins in which markers were assigned

 Written by Marcelo Mollinari

 Escola Superior de Agricultura "Luiz de Queiroz"
 Departamento de Genética - São Paulo, Brazil
 Contact: mmollina@usp.br
 First version: 08/2015
 Last update: 09/2015
 */

#include <algorithm>
#include <iostream>
#include <vector>
#include "find_bins.h"
#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <R_ext/PrtUtil.h>

using namespace std;
using namespace Rcpp;

RcppExport SEXP get_bins(SEXP geno_R, SEXP exact_R)
{
  int exact = Rcpp::as<int>(exact_R);
  Rcpp::NumericMatrix geno = Rcpp::as<Rcpp::NumericMatrix>(geno_R);
  int n_mar = geno.ncol();
  int n_ind = geno.nrow();
  std::vector<int> b_vec(1);
  b_vec[0]=1;
  std::vector<int> b(n_mar);
  std::fill(b.begin(), b.end(), 1);
  int flag = 0, l;
  if(exact)
  {
    for(int i = 0; i < n_mar; i++)
    {
      R_CheckUserInterrupt();
      for(int j = b_vec.size() - 1; j >= 0 ; j--)
      {
        flag=0;
        l=check_occurrence(b, b_vec[j]);
        for(int k = 0; k < n_ind; k++)
        {
          if(geno(k,i)!=geno(k,l))
          {
            flag=1;
            break;
          }
        }
        if(flag==0)
        {
          b[i]=b_vec[j];
          break;
        }
      }
      if(flag==1)
      {
        b[i] = *std::max_element(b_vec.begin(), b_vec.end())+1;
        b_vec.push_back(b[i]);
      }
    }
  }
  else
  {
    for(int i = 0; i < n_mar; i++)
    {
      for(int j = b_vec.size() - 1; j >= 0 ; j--)
      {
        flag=0;
        l=check_occurrence(b, b_vec[j]);
        for(int k = 0; k < n_ind; k++)
        {
          if(geno(k,i)!=geno(k,l) && geno(k,i)!=0 && geno(k,l)!=0)
          {
            flag=1;
            break;
          }
        }
        if(flag==0)
        {
          b[i]=b_vec[j];
          break;
        }
      }
      if(flag==1)
      {
        b[i] = *std::max_element(b_vec.begin(), b_vec.end())+1;
        b_vec.push_back(b[i]);
      }
    }
  }
  return(wrap(b));
}

/*Check if x is contained in v. If yes, it returns its position in v; returns -1 otherwise*/
int check_occurrence(std::vector<int>& v, int x)
{
  for(int i = 0; (unsigned)i < v.size(); i++)
    if(v[i]==x)
      return(i);
        return(-1);
}

