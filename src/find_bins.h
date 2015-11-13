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
 File: find_bins.h
 Description: Header file for 'find_bins' function.
              Search for markers with same inforation and attributes them to bins

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
#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <R_ext/PrtUtil.h>

RcppExport SEXP get_bins(SEXP geno_R, SEXP exact_R, SEXP barWidth_R);
int check_occurrence(std::vector<int>& v, int x);
