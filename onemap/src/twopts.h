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
  File: twopts.h
  Description: Header file for OneMap functions written in C code (two-point analysis)
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

void mdrct2pt(double A[16], double B[16], double res[16]);

void mprod(double *A, int rowA, int colA, double *B, int rowB, int colB, double *res);

void H1(double r, double H[16]);

void H2(double r, double H[16]);

void H3(double r, double H[16]);

void H4(double r, double H[16]);

double log_add(double x, double y);

double log_sub(double x, double y);

void rf_2pt(double *I1, int p1, double *I2, int p2, int *n, int ntot,
              void (*Hcall)(double, double [16]),
			  double D[16], double *rf_assign, double *log_like_assign);

void r2pts(double *I1, int *p1, double *I2, int *p2, int *n, int *ntot,
           double *r, double *log_like, double *posterior, double *LOD);
