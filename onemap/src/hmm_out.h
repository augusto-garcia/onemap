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
  File: hmm_out.h
  Description: Header file for HMM functions written in C code
               Implements the methodology of Hidden Markov Models (HMM)
	       to construct multipoint linkage maps in outcrossing species

  Written by Marcelo Mollinari
  Adapted from hmm_main.h, hmm_f2.h and util.h (found in the R package qtl)
  copyright (c) 2001-10, Karl W Broman                                

  Escola Superior de Agricultura "Luiz de Queiroz"
  Departamento de Genética - São Paulo, Brazil
  Contact: mmollina@usp.br
  First version: 03/02/2009
  Last update: 03/02/2009
*/

double addlog(double a, double b);

void reorg_geno(int n_ind, int n_pos, int *geno, int ***Geno);

void reorg_genoprob(int n_ind, int n_pos, int n_gen, double *genoprob, double ****Genoprob);

void allocate_alpha(int n_pos, int n_gen, double ***alpha);

void allocate_double(int n, double **vector);

void allocate_int(int n, int **vector);

void allocate_dmatrix(int n_row, int n_col, double ***matrix);

void allocate_imatrix(int n_row, int n_col, int ***matrix);

double init_outbred(int true_gen);

double emit_outbred(int obs_gen, int true_gen, double error_prob, int mark_type);

double step_outbred(int gen1, int gen2, int phase, double rf);

double nrec_outbred(int gen1, int gen2, int phase);

void est_map(int n_ind, int n_mar, int *type, int *phase, int n_gen, int *geno, double *rf,
	     double error_prob, double initf(int),
	     double emitf(int, int, double, int),
	     double stepf(int, int, int, double),
	     double nrec(int, int,int),
	     double *loglik, int maxit, double tol,
	     int verbose);

void calc_genoprob(int n_ind, int n_pos, int *type, int *phase, int n_gen, int *geno,
		   double *rf,
		   double error_prob, double *genoprob,
		   double initf(int),
		   double emitf(int, int, double, int),
		   double stepf(int, int, int, double));

void est_map_outbred(int *n_ind, int *n_mar, int *type, int *phase, int *geno, double *rf,
		     double *error_prob, double *loglik, int *maxit,
		     double *tol, int *verbose);

void calc_genoprob_outbred(int *n_ind, int *n_mar, int *type,  int *phase, int *geno,
			   double *rf, double *error_prob, double *genoprob);
