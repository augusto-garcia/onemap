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
 File: vcf2raw.c
 Description: Scan a VCF file and search for markers with potentially known segregation patterns in progeny (given parent samples)
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

/*****************************
 * SEARCH FOR value IN array *
 *****************************/
bool is_val_in_arr(int value, int *array, int arr_length){
  int i;
  for (i = 0; i < arr_length; i++) {
    if (array[i] == value)
      return true;
  }
  return false;
}

/****************************************
 * SEARCH SAMPLE NAMES FOR PARENT NAMES *
 ****************************************/
void get_parents_idx(int n_parent1, int idx_parent1[],
                     int n_parent2, int idx_parent2[],
                     bcf_hdr_t *vcf_hdr, char **parent1, char **parent2) {
  int i, s, n_samples = bcf_hdr_nsamples(vcf_hdr);

  for (i = 0; i < n_parent1; i++) {
    idx_parent1[i] = EMPTY_SAMPLE_INIT;
  }
  for (i = 0; i < n_parent2; i++) {
    idx_parent2[i] = EMPTY_SAMPLE_INIT;
  }

  for (s = 0; s < n_samples; s++) {
    for (i = 0; i < n_parent1; i++) {
      if(!strcmp(vcf_hdr->samples[s], parent1[i])) {
        idx_parent1[i] = s;
        break;
      }
    }
    for (i = 0; i < n_parent2; i++) {
      if(!strcmp(vcf_hdr->samples[s], parent2[i])) {
        idx_parent2[i] = s;
        break;
      }
    }
  }

  // Check for any missing samples
  for (i = 0; i < n_parent1; i++) {
    if(idx_parent1[i] == EMPTY_SAMPLE_INIT) {
      error("Could not find sample %s.\n", parent1[i]);
    }
  }
  for (i = 0; i < n_parent2; i++) {
    if(idx_parent2[i] == EMPTY_SAMPLE_INIT) {
      error("Could not find sample %s.\n", parent2[i]);
    }
  }
}

/******************
 * GET CROSS TYPE *
 ******************/
int get_cross_type(char** cross) {
  int type = cross_NONE;

  if(!strcmp(*cross, "outcross")) {
    type = cross_out;
  }
  else if(!strcmp(*cross, "backcross")) {
    type = cross_bc;
  }
  else if(!strcmp(*cross, "riself") || !strcmp(*cross, "risib")) {
    type = cross_ril;
  }
  else if(!strcmp(*cross, "intercross")) {
    type = cross_f2;
  }

  if (type == cross_NONE) {
    error("Unknown cross type: %s.\n", *cross);
  }

  return type;
}

/**********************************************************
 * CALL PARENT GENOTYPES GIVEN ONE OR MORE SAMPLE INDICES *
 **********************************************************/
void get_consensus_parent_gt(bcf_fmt_t *fmt_ptr, int n_parent, int idx_parent[], int min_class_parent,
                             bool* is_het_parent, bool* is_hom_ref_parent, bool* is_hom_alt_parent) {
  int i, ial, jal;

  // We only consider the following genotypes (https://github.com/samtools/htslib/blob/develop/htslib/vcfutils.h)
  //    0: Homozygous reference allele
  //    1: Homozygous alternative allele (and assume there is only one alternative allele in variant locus)
  //    2: Heterozygous ref/alt
  int count_parent[3] = {0};

  for (i = 0; i < n_parent; i++) {
    int cur_gt = bcf_gt_type(fmt_ptr, idx_parent[i], &ial, &jal);
    if (cur_gt <= 2 && cur_gt >= 0) {
      count_parent[cur_gt]++;
    }
  }

  // First, decide whether parent is heterozygous
  if (count_parent[GT_HET_RA]) {
    if (count_parent[GT_HET_RA] >= min_class_parent) {
      *is_het_parent = true;
    }
  }
  // If not, check whether it is REF allele homozygous or ALT allele homozygous
  else if (count_parent[GT_HOM_RR] && !count_parent[GT_HOM_AA]) {
    if (count_parent[GT_HOM_RR] >= min_class_parent) {
      *is_hom_ref_parent = true;
    }
  }
  else if (!count_parent[GT_HOM_RR] && count_parent[GT_HOM_AA]) {
    if (count_parent[GT_HOM_AA] >= min_class_parent) {
      *is_hom_alt_parent = true;
    }
  }
}

/**********************************************************************************
 * BASED ON CROSS TYPE AND PARENT GENOTYPES, DECIDE ON MARKER SEGREGATION PATTERN *
 **********************************************************************************/
int get_marker_type(char marker_type[MARKER_TYPE_LEN], int cross_type,
                    bool is_het_parent1, bool is_hom_ref_parent1, bool is_hom_alt_parent1,
                    bool is_het_parent2, bool is_hom_ref_parent2, bool is_hom_alt_parent2) {
  int type = marker_NONE;

  switch(cross_type)
  {
  case cross_out:
    if (is_het_parent1) {
      if (is_het_parent2) {
        strcpy(marker_type, "B3.7");
        type = marker_B3;
      }
      else if (is_hom_ref_parent2) {
        strcpy(marker_type, "D1.10");
        type = marker_D_ref;
      }
      else if (is_hom_alt_parent2) {
        strcpy(marker_type, "D1.10");
        type = marker_D_alt;
      }
    }
    else if (is_het_parent2) {
      if (is_hom_ref_parent1) {
        strcpy(marker_type, "D2.15");
        type = marker_D_ref;
      }
      else if (is_hom_alt_parent1) {
        strcpy(marker_type, "D2.15");
        type = marker_D_alt;
      }
    }
    break;
  case cross_bc:
    if (is_hom_ref_parent1 && is_hom_alt_parent2) {
      strcpy(marker_type, "A.H");
      type = marker_BC_ref;
    }
    else if (is_hom_alt_parent1 && is_hom_ref_parent2) {
      strcpy(marker_type, "A.H");
      type = marker_BC_alt;
    }
    break;
  case cross_ril:
    if (is_hom_ref_parent1 && is_hom_alt_parent2) {
      strcpy(marker_type, "A.B");
      type = marker_RI_ref;
    }
    else if (is_hom_alt_parent1 && is_hom_ref_parent2) {
      strcpy(marker_type, "A.B");
      type = marker_RI_alt;
    }
    break;
  case cross_f2:
    if (is_hom_ref_parent1 && is_hom_alt_parent2) {
      // We only consider codominant markers for F2 intercrosses
      strcpy(marker_type, "A.H.B");
      type = marker_F2_ref;
    }
    else if (is_hom_alt_parent1 && is_hom_ref_parent2) {
      strcpy(marker_type, "A.H.B");
      type = marker_F2_alt;
    }
    break;
  default:
    error("Unknown cross type.\n");
  }

  return type;
}

/****************************************************
 * OUTPUT A SINGLE MARKER RECORD TO raw OUTPUT FILE *
 ****************************************************/
void print_record(FILE *temp_f, char *marker_name, char marker_type[MARKER_TYPE_LEN], bcf_fmt_t *fmt_ptr,
                  int n_progeny, int idx_progeny[], const char * const(*type_ptr)[GT_TYPES_LEN]) {
  int i, ial, jal;

  fprintf(temp_f, "*%s\t%s\t", marker_name, marker_type);
  int cur_prog_gt = bcf_gt_type(fmt_ptr, idx_progeny[0], &ial, &jal);
  fprintf(temp_f, "%s", (*type_ptr)[cur_prog_gt]);
  for (i = 1; i < n_progeny; i++) {
    cur_prog_gt = bcf_gt_type(fmt_ptr, idx_progeny[i], &ial, &jal);
    fprintf(temp_f, " %s", (*type_ptr)[cur_prog_gt]);
  }
  fprintf(temp_f, "\n");
}

/**************************
 * PROCESS INPUT VCF FILE *
 **************************/
void vcf2raw(char **filename, char **out_filename, char **cross, int *n_parent1,
             char **parent1, int *n_parent2, char **parent2, double *min_class) {
  // We assume the input file exists (checked in R)
  bcf_sweep_t *in_vcf = bcf_sweep_init(*filename);
  if (in_vcf == NULL) {
    bcf_sweep_destroy(in_vcf);
    error("Could not parse input VCF file.");
  }
  bcf_hdr_t *vcf_hdr = bcf_sweep_hdr(in_vcf);

  // Get reference sequence IDs
  int n_seq = 0;
  const char **seq_names = NULL;
  seq_names = bcf_hdr_seqnames(vcf_hdr, &n_seq);
  if (seq_names == NULL || n_seq == 0) {
    free(seq_names);
    error("Could not correctly parse sequence names in VCF file. Is the input file tabix indexed?\n");
  }

  // Map parent names to sample indices
  int idx_parent1[*n_parent1];
  int idx_parent2[*n_parent2];
  get_parents_idx(*n_parent1, idx_parent1, *n_parent2, idx_parent2, vcf_hdr, parent1, parent2);

  // Get progeny sample indices (all samples that are not set as parents)
  int n_samples = bcf_hdr_nsamples(vcf_hdr);
  int n_progeny = n_samples - *n_parent1 - *n_parent2;
  if (n_progeny == 0) {
    error("Input file must contain at least one progeny individual.");
  }
  int idx_progeny[n_progeny];
  int i = 0, s;
  for (s = 0; s < n_samples; s++) {
    if (!is_val_in_arr(s, idx_parent1, *n_parent1)) {
      if (!is_val_in_arr(s, idx_parent2, *n_parent2)) {
        idx_progeny[i++] = s;
      }
    }
  }

  // Minimum count to assign parent genotype
  int min_class_parent1 = (int)ceil(*min_class * *n_parent1);
  int min_class_parent2 = (int)ceil(*min_class * *n_parent2);

  // Convert cross type
  int cross_type = get_cross_type(cross);

  // We need to write to a temporary file, because the number of markers in the header is unknown
  FILE *temp_f;
  char temp_filename[] = "tmp_raw_XXXXXX";
  int temp_fd;
  temp_fd = mkstemp(temp_filename);
  if (temp_fd == -1) {
    error("Could not open temporary output file.\n");
  }
  unlink(temp_filename);
  temp_f = fdopen(temp_fd, "w+");
  if (temp_f == NULL) {
    error("Could not open temporary output file.\n");
  }

  // CHROM and POS fields will be placed at the end of the output file
  int marker_count = 0;
  int * chrom = malloc(MAX_VARIANTS * sizeof(int));
  if (chrom == NULL) {
    error("Could not allocate vector.\n");
  }
  int * pos = malloc(MAX_VARIANTS * sizeof(int));
  if (pos == NULL) {
    error("Could not allocate vector.\n");
  }

  // Mapping of VCF genotypes to ONEMAP genotypes
  const char * const D_BC_ref[GT_TYPES_LEN] = { "a", "-", "ab", "-", "-", "-", "-" };
  const char * const D_BC_alt[GT_TYPES_LEN] = { "-", "a", "ab", "-", "-", "-", "-" };
  const char * const RI_ref[GT_TYPES_LEN] = { "a", "b", "-", "-", "-", "-", "-" };
  const char * const RI_alt[GT_TYPES_LEN] = { "b", "a", "-", "-", "-", "-", "-" };
  const char * const B3_F2_ref[GT_TYPES_LEN] = { "a", "b", "ab", "-", "-", "-", "-" };
  const char * const B3_F2_alt[GT_TYPES_LEN] = { "b", "a", "ab", "-", "-", "-", "-" };

  // Scan all records in VCF file and print valid markers to output
  bcf1_t *record;
  int32_t *GTs = NULL;
  int nGT_arr = 0;

  while ((record = bcf_sweep_fwd(in_vcf)) && marker_count < MAX_VARIANTS) {
    // We only consider biallelic SNP and INDEL markers
    int var_type = bcf_get_variant_types(record);
    if ((var_type == VCF_SNP || var_type == VCF_INDEL) && record->n_allele == 2) {
      int nGTs = bcf_get_format_int32(vcf_hdr, record, "GT", &GTs, &nGT_arr);
      // We only consider diploid variants (number of alleles in genotypes == 2)
      nGTs /= n_samples;
      if (nGTs == 2) {

        bcf_fmt_t *fmt_ptr = bcf_get_fmt(vcf_hdr, record, "GT");

        // First, check which parents are heterozygous or homozygous (REF or ALT allele)
        bool is_het_parent1 = false, is_hom_ref_parent1 = false, is_hom_alt_parent1 = false;
        get_consensus_parent_gt(fmt_ptr, *n_parent1, idx_parent1, min_class_parent1, &is_het_parent1,
                                &is_hom_ref_parent1, &is_hom_alt_parent1);
        bool is_het_parent2 = false, is_hom_ref_parent2 = false, is_hom_alt_parent2 = false;
        get_consensus_parent_gt(fmt_ptr, *n_parent2, idx_parent2, min_class_parent2, &is_het_parent2,
                                &is_hom_ref_parent2, &is_hom_alt_parent2);

        // Convert to appropriate marker type
        char marker_type[MARKER_TYPE_LEN];
        int type = get_marker_type(marker_type, cross_type,
                                   is_het_parent1, is_hom_ref_parent1, is_hom_alt_parent1,
                                   is_het_parent2, is_hom_ref_parent2, is_hom_alt_parent2);

        const char * const(*type_ptr)[GT_TYPES_LEN];
        bool valid_marker = true;
        switch(type)
        {
        case marker_B3:
        case marker_F2_ref:
          type_ptr = &B3_F2_ref;
          break;
        case marker_F2_alt:
          type_ptr = &B3_F2_alt;
          break;
        case marker_D_ref:
        case marker_BC_ref:
          type_ptr = &D_BC_ref;
          break;
        case marker_D_alt:
        case marker_BC_alt:
          type_ptr = &D_BC_alt;
          break;
        case marker_RI_ref:
          type_ptr = &RI_ref;
          break;
        case marker_RI_alt:
          type_ptr = &RI_alt;
          break;
        default:
          valid_marker = false;
        }

        if (valid_marker) {
          // Store CHROM and POS fields for valid markers
          chrom[marker_count] = record->rid;
          pos[marker_count] = record->pos + 1;

          // Check if marker name exists; if negative, create one
          char *marker_name = record->d.id;
          if (!strcmp(marker_name, ".")) {
            sprintf(marker_name, "%s.%d", seq_names[chrom[marker_count]], pos[marker_count]);
          }

          // Output variant in ONEMAP format to temporary file
          print_record(temp_f, marker_name, marker_type, fmt_ptr, n_progeny, idx_progeny, type_ptr);

          marker_count++;
        }
      }
    }
  }

  // Write final output file header
  FILE *final_f = fopen(*out_filename, "w");
  if (final_f == NULL) {
    error("Could not open output file.\n");
  }
  fprintf(final_f, "data type %s\n", *cross);
  // The next header line contains the following information: number of individuals, number of markers, 1 for the presence of CHROM information, 1 for the presence of POS information and 0 for the absence of phenotypes (these need to be manually included later)
  fprintf(final_f, "%d %d 1 1 0\n", n_progeny, marker_count);

  // Copy marker data from temporary file to final file
  rewind(temp_f);
  char buf[BUFSIZ];
  size_t size;
  while ((size = fread(buf, 1, BUFSIZ, temp_f))) {
    fwrite(buf, 1, size, final_f);
  }

  // Write CHROM and POS data to output file
  if (marker_count) {
    fprintf(final_f, "*CHROM\t");
    fprintf(final_f, "%s", seq_names[chrom[0]]);
    for (i = 1; i < marker_count; i++) {
      fprintf(final_f, " %s", seq_names[chrom[i]]);
    }
    fprintf(final_f, "\n*POS\t");
    fprintf(final_f, "%d", pos[0]);
    for (i = 1; i < marker_count; i++) {
      fprintf(final_f, " %d", pos[i]);
    }
  }

  // Clean-up
  free(chrom);
  free(pos);

  free(GTs);
  bcf_sweep_destroy(in_vcf);

  fclose(temp_f);
  close(temp_fd);
  fclose(final_f);
}
