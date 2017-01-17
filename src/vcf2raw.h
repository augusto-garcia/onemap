/*
  OneMap: software for genetic mapping in outcrossing species
  Copyright (C) 2007-2017 Gabriel R A Margarido and Marcelo Mollinari

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
  File: vcf2raw.h
  Description: Header file for 'vcf2raw' functions
               Scan a VCF file and search for markers with potentially known segregation patterns in progeny (given parent samples)
               Output a raw file with valid markers in OneMap format

  Written by Gabriel R A Margarido

  Escola Superior de Agricultura "Luiz de Queiroz"
  Departamento de Genética - São Paulo, Brazil
  Contact: gramarga@usp.br
  First version: 10/2015
  Last update: 12/2015
*/

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <vcf_sweep.h>
#include <vcfutils.h>
#include <R.h>

#define EMPTY_SAMPLE_INIT -1
#define POS_LEN 15
#define MARKER_TYPE_LEN 6
#define GT_TYPES_LEN 7
#define MAX_VARIANTS 100000000

#define cross_NONE 0
#define cross_out 1
#define cross_bc 2
#define cross_ril 3
#define cross_f2 4

#define marker_NONE 0
#define marker_B3 1
#define marker_D_alt 2
#define marker_D_ref 3
#define marker_BC_ref 4
#define marker_BC_alt 5
#define marker_RI_ref 6
#define marker_RI_alt 7
#define marker_F2_ref 8
#define marker_F2_alt 9

bool is_val_in_arr(int value, int *array, int arr_length);

void get_parents_idx(int n_parent1, int idx_parent1[],
		     int n_parent2, int idx_parent2[],
		     bcf_hdr_t *vcf_hdr, char **parent1, char **parent2);

int get_cross_type(char** cross);

void get_consensus_parent_gt(bcf_fmt_t *fmt_ptr, int n_parent, int idx_parent[], int min_class_parent,
			     bool* is_het_parent, bool* is_hom_ref_parent, bool* is_hom_alt_parent);

int get_marker_type(char marker_type[MARKER_TYPE_LEN], int cross_type,
                    bool is_het_parent1, bool is_hom_ref_parent1, bool is_hom_alt_parent1,
                    bool is_het_parent2, bool is_hom_ref_parent2, bool is_hom_alt_parent2);

void print_record(FILE *temp_f, char *marker_name, char marker_type[MARKER_TYPE_LEN], bcf_fmt_t *fmt_ptr,
		  int n_progeny, int idx_progeny[], const char * const(*type_ptr)[GT_TYPES_LEN]);

void vcf2raw(char **filename, char **out_filename, char **cross, int *n_parent1,
             char **parent1, int *n_parent2, char **parent2, double *min_class);
