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
  File: threepts.h
  Description: Header file for OneMap functions written in C code (three-point analysis)
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

void mdrct3pt(double A[64], double B[64], double res[64]);

void mkron(double *A, int rowA, int colA, double *B, int rowB, int colB, double *res);

void H11(double g00, double g01, double g10, double g11, double H[64]);

void H12(double g00, double g01, double g10, double g11, double H[64]);

void H13(double g00, double g01, double g10, double g11, double H[64]);

void H14(double g00, double g01, double g10, double g11, double H[64]);

void H21(double g00, double g01, double g10, double g11, double H[64]);

void H22(double g00, double g01, double g10, double g11, double H[64]);

void H23(double g00, double g01, double g10, double g11, double H[64]);

void H24(double g00, double g01, double g10, double g11, double H[64]);

void H31(double g00, double g01, double g10, double g11, double H[64]);

void H32(double g00, double g01, double g10, double g11, double H[64]);

void H33(double g00, double g01, double g10, double g11, double H[64]);

void H34(double g00, double g01, double g10, double g11, double H[64]);

void H41(double g00, double g01, double g10, double g11, double H[64]);

void H42(double g00, double g01, double g10, double g11, double H[64]);

void H43(double g00, double g01, double g10, double g11, double H[64]);

void H44(double g00, double g01, double g10, double g11, double H[64]);

void rf_3pt(double *I1, int p1, double *I2, int p2, double *I3, int p3, int *n, int ntot, void (*Hcall)(double, double, double, double, double [64]), double G00[64], double G01[64], double G10[64], double G11[64], double *theta12_assign, double *theta23_assign, double *theta13_assign, double *log_like_assign);

void r3pts(double *I1, int *p1, double *I2, int *p2, double *I3, int *p3, int *n, int *ntot, double *theta12, double *theta23, double *theta13, double *log_like, double *posterior, double *LOD);
