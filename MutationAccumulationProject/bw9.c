// bw9.c
// Authors: Ruchira S. Datta, Trevor Graham
// Copyright (c) 2012, Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// o Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// o Redistributions in binary form mus reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// o Neither the name of the University of California, San Francisco nor the
// names of its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS" AND ANY 
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PUPROSE ARE
// DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMTIED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

/*
 
 bw2: included support for readindata
 bw3: switched to multinomial, corrected mutation mistake (cells could get mutated multiple times)
 bw4: reverted to type int, corrected mutation mistake (cells could get mutated multiple times)
 bw5: moved to unsigned int, corrected mutation mistakes, allowed multiple mutations in a single cell
 bw6: corrected and simplified subclone
 bw7: multiple mutations per clone per generation
 bw8: added thresholding for cancer: cancer cells need to constitute X proportion of the tumor
 bw9: added passenger mutation class;
      used SPRNG 2.0 with MPI for parallel random number generation; 
      used gsl vector, gsl matrix, and marray;
      used forward simulation;
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#define SIMPLE_SPRNG
#define USE_MPI
#include <gsl/gsl_rng.h>
#include "gsl-sprng.h"
#include <gsl/gsl_randist.h>
#include "readindata.h"
#include <gsl/gsl_vector_uint.h>
#include <gsl/gsl_matrix.h>
#include "marray_uint.h"
#include <gsl/gsl_math.h>

#define MAXMUTSPERCELL 4

#define CLASSBLOCKSIZE  256

struct mutant_class_struct {
  size_t indices[4];
};

typedef struct mutant_class_struct mutant_class_indices;

struct diversity_indices_struct {
  int richness;
  double shannon_index;
  double simpson_index;
};

typedef struct diversity_indices_struct DiversityIndices;

//#define DEBUG

void compute_fitness_coeff(gsl_matrix *F, int Md, int Mh);
void fill_powers(double *powers, double alpha, int highest_power);
void fill_prob_of_cells_acquiring_num_muts(
                                  double **prob_of_cells_acquiring_num_muts,
                                  int total_num_loci,
                                  double mu,
                                  double *mupower,
                                  double *oneminusmupower,
                                  double *small_reciprocal_factorial);
void fill_reciprocal_lookup(double *reciprocal_num_remaining_loci, int
                            total_num_loci);

unsigned int do_cycle(gsl_rng *gslr, marray_uint *counts, 
  gsl_matrix *F, 
  unsigned int N, int Md, int Mh, int Mp, int Mm,
  int *min_occupied_D_ptr, int *min_occupied_H_ptr,
  int *min_occupied_P_ptr, int *min_occupied_M_ptr,
  int *max_occupied_D_ptr, int *max_occupied_H_ptr,
  int *max_occupied_P_ptr, int *max_occupied_M_ptr,
  int *num_potentially_occupied_mutant_classes_ptr,
  mutant_class_indices **potentially_occupied_mutant_classes_ptr,
  marray_uint *index_of_potentially_occupied_mutant_class,
  double **mutant_class_probabilities_ptr,
  unsigned int **new_mutant_class_counts_ptr,
  marray_uint *delta_counts,
  double **prob_of_cells_acquiring_num_muts_lo, 
  double **prob_of_cells_acquiring_num_muts_hi, 
  unsigned int num_cells_acquiring_num_muts[], 
  double prob_that_mutation_falls_in_class[],
  double *small_reciprocal_factorial,
  double *reciprocal_num_remaining_loci,
  size_t *indices, size_t *delta_indices, size_t *new_indices);
unsigned int grow_total_counts(gsl_rng *gslr, marray_uint *counts, gsl_matrix *F, int Md, int Mh, int Mp, int Mm, size_t *indices);
void mutate_counts(gsl_rng *gslr, marray_uint *counts, int Md, int Mh, int Mp, int Mm, 
    int *max_occupied_D_ptr, int *max_occupied_H_ptr,
    int *max_occupied_P_ptr, int *max_occupied_M_ptr,
    marray_uint *delta_counts,
    double **prob_of_cells_acquiring_num_muts_lo,
    double **prob_of_cells_acquiring_num_muts_hi,
    unsigned int num_cells_acquiring_num_muts[],
    double prob_that_mutation_falls_in_class[], 
    size_t *indices, size_t *delta_indices, size_t *new_indices);
void mutate_counts_in_class(gsl_rng *gslr, marray_uint *counts, int Md, int Mh, int Mp, int Mm, 
    int *max_occupied_D_ptr, int *max_occupied_H_ptr,
    int *max_occupied_P_ptr, int *max_occupied_M_ptr,
    marray_uint *delta_counts,
    double **prob_of_cells_acquiring_num_muts,
    unsigned int num_cells_acquiring_num_muts[],
    double prob_that_mutation_falls_in_class[], 
    size_t *indices, size_t *delta_indices, size_t *new_indices, 
    int mutation_rate_level);
void distribute_mutations(gsl_rng *gslr, marray_uint *delta_counts, size_t *delta_indices,
                          double *prob_that_mutation_falls_in_class, 
                          unsigned int num_cells_acquiring_num_mutations, 
                          unsigned int num_mutations_left);
void update_class_counts(gsl_rng *gslr, int N, marray_uint *counts, 
                        gsl_matrix *F,
                        int Md, int Mh, int Mp, int Mm,
                        int *min_occupied_D_ptr, int *min_occupied_H_ptr,
                        int *min_occupied_P_ptr, int *min_occupied_M_ptr,
                        int *max_occupied_D_ptr, int *max_occupied_H_ptr,
                        int *max_occupied_P_ptr, int *max_occupied_M_ptr,
                        int *num_potentially_occupied_mutant_classes_ptr,
                        mutant_class_indices **potentially_occupied_mutant_classes_ptr,
                        marray_uint *index_of_potentially_occupied_mutant_class,
                        double **mutant_class_probabilities_ptr,
                        unsigned int **new_mutant_class_counts_ptr,
                        double **prob_of_cells_acquiring_num_muts_lo,
                        double **prob_of_cells_acquiring_num_muts_hi,
                        double *small_reciprocal_factorial,
                        double *reciprocal_num_remaining_loci,
                        size_t *indices, size_t *delta_indices, 
                        size_t *new_indices);

int compute_most_aggressive_mutant_class(marray_uint *counts, 
                                        int min_occupied_D, int min_occupied_H,
                                        int min_occupied_P, int min_occupied_M,
                                        int max_occupied_D, int max_occupied_H,
                                        int max_occupied_P, int max_occupied_M,
                                        int *class, size_t *indices);
int compute_most_mutator_mutant_class(marray_uint *counts, 
                                      int min_occupied_D, int min_occupied_H,
                                      int min_occupied_P, int min_occupied_M,
                                      int max_occupied_D, int max_occupied_H,
                                      int max_occupied_P, int max_occupied_M,
                                      int *class, size_t *indices);
int compute_least_domestic_mutant_class(marray_uint *counts, 
                                      int min_occupied_D, int min_occupied_H, 
                                      int min_occupied_P, int min_occupied_M,
                                      int max_occupied_D, int max_occupied_H, 
                                      int max_occupied_P, int max_occupied_M,
                                      int *class, size_t *indices);
int compute_most_hitchhiking_mutant_class(marray_uint *counts, 
                                      int min_occupied_D, int min_occupied_H, 
                                      int min_occupied_P, int min_occupied_M,
                                      int max_occupied_D, int max_occupied_H, 
                                      int max_occupied_P, int max_occupied_M,
                                      int *class, size_t *indices);
int largest_mutant_class(marray_uint *counts, 
                        int min_occupied_D, int min_occupied_H,
                        int min_occupied_P, int min_occupied_M,
                        int max_occupied_D, int max_occupied_H,
                        int max_occupied_P, int max_occupied_M,
                        int *class, size_t *indices);
int compute_number_of_cells_classed_as_cancerous(marray_uint *counts, 
                                                int min_occupied_D, int min_occupied_H,
                                                int min_occupied_P, int min_occupied_M,
                                                int max_occupied_D, int max_occupied_H,
                                                int max_occupied_P, int max_occupied_M,
                                                size_t *indices);
int compute_number_of_mutator_cells(marray_uint *counts, 
                                    int min_occupied_D, int min_occupied_H,
                                    int min_occupied_P, int min_occupied_M,
                                    int max_occupied_D, int max_occupied_H,
                                    int max_occupied_P, int max_occupied_M,
                                    size_t *indices);
void compute_marginalized_driver_counts(marray_uint *counts, int Md, 
                                    int min_occupied_D, int min_occupied_H,
                                    int min_occupied_P, int min_occupied_M,
                                    int max_occupied_D, int max_occupied_H,
                                    int max_occupied_P, int max_occupied_M,
                                    size_t *indices,
                                    marray_uint *marginalized_driver_counts, 
                                    size_t *marginalized_indices, FILE *fp);
void compute_marginalized_housekeeper_counts(marray_uint *counts, int Mh, 
                                    int min_occupied_D, int min_occupied_H,
                                    int min_occupied_P, int min_occupied_M,
                                    int max_occupied_D, int max_occupied_H,
                                    int max_occupied_P, int max_occupied_M,
                                    size_t *indices,
                                    marray_uint *marginalized_housekeeper_counts, 
                                    size_t *marginalized_indices, FILE *fp);
void compute_marginalized_passenger_counts(marray_uint *counts, int Mp, 
                                    int min_occupied_D, int min_occupied_H,
                                    int min_occupied_P, int min_occupied_M,
                                    int max_occupied_D, int max_occupied_H,
                                    int max_occupied_P, int max_occupied_M,
                                    size_t *indices,
                                    marray_uint *marginalized_passenger_counts, 
                                    size_t *marginalized_indices, FILE *fp);
void compute_marginalized_mutator_counts(marray_uint *counts, int Mm, 
                                    int min_occupied_D, int min_occupied_H,
                                    int min_occupied_P, int min_occupied_M,
                                    int max_occupied_D, int max_occupied_H,
                                    int max_occupied_P, int max_occupied_M,
                                    size_t *indices,
                                    marray_uint *marginalized_mutator_counts, 
                                    size_t *marginalized_indices, FILE *fp);
void compute_averages(marray_uint *counts, unsigned int N, 
                      int min_occupied_D, int min_occupied_H,
                      int min_occupied_P, int min_occupied_M,
                      int max_occupied_D, int max_occupied_H,
                      int max_occupied_P, int max_occupied_M,
                      gsl_matrix *dataM, int index, size_t *indices);
void report_data(int process_number, marray_uint *counts, unsigned int N, 
                  int rep, int gen, 
                  int min_occupied_D, int min_occupied_H,
                  int min_occupied_P, int min_occupied_M,
                  int max_occupied_D, int max_occupied_H,
                  int max_occupied_P, int max_occupied_M,
                  int time_first_sig_mut, 
                  int time_first_sig_ca, size_t *indices, 
                  DiversityIndices *diversity_indices_ptr, FILE *fp);
void compute_diversity_indices(marray_uint *counts,
                              int min_occupied_D, int min_occupied_H, 
                              int min_occupied_P, int min_occupied_M,
                              int max_occupied_D, int max_occupied_H, 
                              int max_occupied_P, int max_occupied_M,
                              size_t *indices,
                              double reciprocal_population_size,
                              DiversityIndices *diversity_indices_ptr);

inline int max(int a, int b, int c, int d) {
  return GSL_MAX_INT(GSL_MAX_INT(GSL_MAX_INT(a, b), c), d);
}

int main(int argc, char *argv[]) {
    
  // random number generator
  gsl_rng *gslr;

  // read in initialConditions
  readics(argc, argv);
    
  gsl_matrix *F;
  marray_uint *counts;
  marray_uint *delta_counts;
	unsigned int N;
	int Md,Mh,Mp,Mm,cancerP,mutP;
	int time_first_sig_mut = ics.gens + 1;
	int time_first_sig_ca = ics.gens + 1;
	double beta,prop;
  size_t tensor_dimensions[4];
  size_t indices[4];
  size_t delta_indices[4];
  size_t new_indices[4];

  int process_number;
    
	int class[4];
  int *intptr;

	int i,j,rep;
  double denom, reciprocal_denom;
  int total_num_loci;
  double **prob_of_cells_acquiring_num_muts_lo;
  double **prob_of_cells_acquiring_num_muts_hi;
  unsigned int num_cells_acquiring_num_muts[MAXMUTSPERCELL+1];
  double prob_that_mutation_falls_in_class[4];
  double *reciprocal_num_remaining_loci;

  double mupower[MAXMUTSPERCELL+1];
  double *oneminusmupower;
  int factorial;
  double small_reciprocal_factorial[MAXMUTSPERCELL+1];
  double *small_reciprocal_factorial_ptr;

  gsl_vector_uint *n;
	gsl_matrix *dataM;
	gsl_vector_uint *mutn;

  marray_uint *marginalized_driver_counts_per_generation;
  marray_uint *marginalized_housekeeper_counts_per_generation;
  marray_uint *marginalized_passenger_counts_per_generation;
  marray_uint *marginalized_mutator_counts_per_generation;
  size_t marginalized_dimensions[2];
  size_t marginalized_indices[2];
  int min_occupied_D, min_occupied_H, min_occupied_P, min_occupied_M;
  int max_occupied_D, max_occupied_H, max_occupied_P, max_occupied_M;

  int num_potentially_occupied_mutant_classes;
  mutant_class_indices *potentially_occupied_mutant_classes;
  double *mutant_class_probabilities;
  unsigned int *new_mutant_class_counts;
  marray_uint *index_of_potentially_occupied_mutant_class;
  DiversityIndices diversity_indices;

  // Having less than MAXMUTSPERCELL loci causes a bug, so force there to be at
  // least this many
  ics.Md = GSL_MAX_INT(ics.Md, MAXMUTSPERCELL+1);
  ics.Mh = GSL_MAX_INT(ics.Mh, MAXMUTSPERCELL+1);
  ics.Mp = GSL_MAX_INT(ics.Mp, MAXMUTSPERCELL+1);
  ics.Mm = GSL_MAX_INT(ics.Mm, MAXMUTSPERCELL+1);
	Md=ics.Md;
	Mh=ics.Mh;
  Mp=ics.Mp;
	Mm=ics.Mm;
	beta=ics.beta;

  small_reciprocal_factorial_ptr = small_reciprocal_factorial;
  *small_reciprocal_factorial_ptr = 1.0;
  for (i = 1, factorial = 1, ++small_reciprocal_factorial_ptr; 
      i <= MAXMUTSPERCELL; 
      ++i, factorial *= i, ++small_reciprocal_factorial_ptr) {
    *small_reciprocal_factorial_ptr = 1.0 / factorial;
  }

  total_num_loci = (Md + 1) + (Mh + 1) + (Mp + 1) + (Mm + 1);
  prob_of_cells_acquiring_num_muts_lo 
    = (double **)malloc( (total_num_loci + 1) * sizeof(double *));
  if (!prob_of_cells_acquiring_num_muts_lo) {
    printf("Unable to allocate memory for prob_of_cells_acquiring_num_muts_lo\n");
    exit(1);
  }
  for (i = 0; i <= total_num_loci; ++i) {
    prob_of_cells_acquiring_num_muts_lo[i] 
      = (double *)malloc((MAXMUTSPERCELL + 1) * sizeof(double));
    if (!prob_of_cells_acquiring_num_muts_lo[i]) {
      printf("Unable to allocate memory for prob_of_cells_acquiring_num_muts_lo[%d]\n", i);
      exit(1);
    }
  }
  prob_of_cells_acquiring_num_muts_hi 
    = (double **)malloc( (total_num_loci + 1) * sizeof(double *));
  if (!prob_of_cells_acquiring_num_muts_hi) {
    printf("Unable to allocate memory for prob_of_cells_acquiring_num_muts_hi\n");
    exit(1);
  }
  for (i = 0; i <= total_num_loci; ++i) {
    prob_of_cells_acquiring_num_muts_hi[i] 
      = (double *)malloc((MAXMUTSPERCELL + 1) * sizeof(double));
    if (!prob_of_cells_acquiring_num_muts_hi[i]) {
      printf("Unable to allocate memory for prob_of_cells_acquiring_num_muts_hi[%d]\n", i);
      exit(1);
    }
  }
  oneminusmupower = (double *)malloc((total_num_loci + 1) * sizeof(double));
  if (!oneminusmupower) {
    printf("Unable to allocate memory for oneminusmupower\n");
    exit(1);
  }
  fill_prob_of_cells_acquiring_num_muts(prob_of_cells_acquiring_num_muts_lo,
                                        total_num_loci,
                                        ics.muB,
                                        mupower,
                                        oneminusmupower,
                                        small_reciprocal_factorial);
  fill_prob_of_cells_acquiring_num_muts(prob_of_cells_acquiring_num_muts_hi,
                                        total_num_loci,
                                        ics.sM,
                                        mupower,
                                        oneminusmupower,
                                        small_reciprocal_factorial);
  reciprocal_num_remaining_loci 
    = (double *)malloc((total_num_loci + 1) * sizeof(double));
  if (!reciprocal_num_remaining_loci) {
    printf("Unable to allocate memory for reciprocal_num_remaining_loci\n");
    exit(1);
  }
  fill_reciprocal_lookup(reciprocal_num_remaining_loci, total_num_loci);

  prob_that_mutation_falls_in_class[0] = ics.probD;
  prob_that_mutation_falls_in_class[1] = ics.probH;
  prob_that_mutation_falls_in_class[2] 
    = (1.0 - ics.probD - ics.probH - ics.probM);
  prob_that_mutation_falls_in_class[3] = ics.probM;
	
  printf("#%d Initializing MPI...", process_number);
  // Initialize MPI for the SPRNG 2.0 scalable parallel random number
  // generator
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_number);
  printf("process %d, done\n", process_number);
  if (ics.run_control_file[0] != '\0') {
    printf("#process %d using line %d from %s\n",
      process_number, 
      ics.run_number, ics.run_control_file);
  }
  // Only print the parameters to the output file once.
  if (process_number == 0) {
    printParameters();
  }

  gslr = gsl_rng_alloc (gsl_rng_sprng20);
    
	// assign memory
  tensor_dimensions[0] = (size_t)(Md+1);
  tensor_dimensions[1] = (size_t)(Mh+1);
  tensor_dimensions[2] = (size_t)(Mp+1);
  tensor_dimensions[3] = (size_t)(Mm+1);
	counts = marray_uint_calloc(4, tensor_dimensions);
  index_of_potentially_occupied_mutant_class = marray_uint_calloc(4,
                                                          tensor_dimensions);
  tensor_dimensions[0] = (size_t)(MAXMUTSPERCELL+1);
  tensor_dimensions[1] = (size_t)(MAXMUTSPERCELL+1);
  tensor_dimensions[2] = (size_t)(MAXMUTSPERCELL+1);
  tensor_dimensions[3] = (size_t)(MAXMUTSPERCELL+1);
  delta_counts = marray_uint_calloc(4, tensor_dimensions);

	F = gsl_matrix_calloc(Md+1,Mh+1);
	
	n = gsl_vector_uint_calloc((Md+1)*(Mh+1)*(Mp+1)*(Mm+1)); // for subclone growth 
	
	mutn = gsl_vector_uint_calloc(ics.mut_max);
	dataM = gsl_matrix_calloc(ics.gens/100+1,10);

  marginalized_dimensions[0] = (size_t)ics.gens;
  marginalized_dimensions[1] = (size_t)(Md+1);
  marginalized_driver_counts_per_generation 
    = marray_uint_calloc(2, marginalized_dimensions);
  marginalized_dimensions[1] = (size_t)(Mh+1);
  marginalized_housekeeper_counts_per_generation
    = marray_uint_calloc(2, marginalized_dimensions);
  marginalized_dimensions[1] = (size_t)(Mp+1);
  marginalized_passenger_counts_per_generation
    = marray_uint_calloc(2, marginalized_dimensions);
  marginalized_dimensions[1] = (size_t)(Mm+1);
  marginalized_mutator_counts_per_generation
    = marray_uint_calloc(2, marginalized_dimensions);
	
  num_potentially_occupied_mutant_classes = (MAXMUTSPERCELL + 1)^4;
  potentially_occupied_mutant_classes = (mutant_class_indices *)malloc(
        num_potentially_occupied_mutant_classes * sizeof(mutant_class_indices));
  if (!potentially_occupied_mutant_classes) {
    printf("Unable to allocate memory for potentially_occupied_mutant_classes\n");
    exit(1);
  }
  mutant_class_probabilities = (double *)malloc(
            num_potentially_occupied_mutant_classes * sizeof(double));
  if (!mutant_class_probabilities) {
    printf("Unable to allocate memory for mutant_class_probabilities\n");
    exit(1);
  }
  new_mutant_class_counts = (unsigned int *)malloc(
            num_potentially_occupied_mutant_classes * sizeof(unsigned int));
  if (!new_mutant_class_counts) {
    printf("Unable to allocate memory for new_mutant_class_counts\n");
    exit(1);
  }

	// initialize
	// calculate fitness coefficients
	compute_fitness_coeff(F,Md,Mh);

  // initialize class counts to zero
  for (i = 0, intptr = class; i < 4; ++i) {
    *intptr++ = 0;
  }
	
  char filename[255];
  sprintf(filename, "%s_generational_proc_%05d", ics.filename, process_number);
	FILE *fp;
	
	if (!(fp = fopen(filename, "w"))) {
		printf("#*** error! couldn't open %s for writing\n",ics.filename);
		exit(0);
	}
  FILE *marginalized_driver_count_fp = NULL;
  FILE *marginalized_housekeeper_count_fp = NULL;
  FILE *marginalized_passenger_count_fp = NULL;
  FILE *marginalized_mutator_count_fp = NULL;
  if (ics.output_marginalized_counts==1) {
    sprintf(filename, "%s-marginalized-driver-counts_proc_%05d", ics.filename,
            process_number);
    if (!(marginalized_driver_count_fp = fopen(filename, "w"))) {
      printf("#*** error! couldn't open %s for writing\n",filename);
      exit(0);
    }
    sprintf(filename, "%s-marginalized-housekeeper-counts_proc_%05d", 
          ics.filename, process_number);
    if (!(marginalized_housekeeper_count_fp = fopen(filename, "w"))) {
      printf("#*** error! couldn't open %s for writing\n",filename);
      exit(0);
    }
    sprintf(filename, "%s-marginalized-passenger-counts_proc_%05d", 
            ics.filename, process_number);
    if (!(marginalized_passenger_count_fp = fopen(filename, "w"))) {
      printf("#*** error! couldn't open %s for writing\n",filename);
      exit(0);
    }
    sprintf(filename, "%s-marginalized-mutator-counts_proc_%05d", 
            ics.filename, process_number);
    if (!(marginalized_mutator_count_fp = fopen(filename, "w"))) {
      printf("#*** error! couldn't open %s for writing\n",filename);
      exit(0);
    }
  }
  FILE *average_fp;
  sprintf(filename, "%s-averages_proc_%05d", ics.filename, process_number);
  if (!(average_fp = fopen(filename, "w"))) {
    printf("#*** error! couldn't open %s for writing\n",filename);
    exit(0);
  }
  FILE *got_cancer_fp;
  sprintf(filename, "%s-got_cancer_proc_%05d", ics.filename, process_number);
  if (!(got_cancer_fp = fopen(filename, "w"))) {
    printf("#*** error! couldn't open %s for writing\n",filename);
    exit(0);
  }

	printParametersToFile(fp);
  fprintf(got_cancer_fp,"GotCancerP,");
  FILE *fileptr;
	// print header to files
  for (i = 0; i < 2; ++i) {
    if (i == 0) {
      fileptr = fp;
    } else {
      fileptr = got_cancer_fp;
    }
    fprintf(fileptr,"proc_num,rc_line,");
    fprintf(fileptr,"maxgens,rep,gen,N,");
    fprintf(fileptr,"mD,mH,mP,mM,");
    fprintf(fileptr,"min_occupied_D,max_occupied_D,");
    fprintf(fileptr,"min_occupied_H,max_occupied_H,");
    fprintf(fileptr,"min_occupied_P,max_occupied_P,");
    fprintf(fileptr,"min_occupied_M,max_occupied_M,");
    fprintf(fileptr,"sD,sH,sM,muB,cancer_muts,");
    fprintf(fileptr,"mean_mD,mean_mH,mean_mP,mean_mM,");
    fprintf(fileptr,"num_in_largest_class,largest_class_mD,largest_class_mH,largest_class_mP,largest_class_mM,");
    fprintf(fileptr,"num_in_most_aggressive_class,most_aggressive_class_mD,most_aggressive_class_mH,most_aggressive_class_mP,most_aggressive_class_mM,");
    fprintf(fileptr,"num_in_most_mutator_class,most_mutator_class_mD,most_mutator_class_mH,most_mutator_class_mP,most_mutator_class_mM,");
    fprintf(fileptr,"num_in_most_housekeeper_class,most_housekeeper_class_mD,most_housekeeper_class_mH,most_housekeeper_class_mP,most_housekeeper_class_mM,");
    fprintf(fileptr,"num_in_most_passenger_class,most_passenger_class_mD,most_passenger_class_mH,most_passenger_class_mP,most_passenger_class_mM,");
    fprintf(fileptr,"time_first_sig_mut,prop_final_mut,time_first_sig_ca,prop_final_can,");
    fprintf(fileptr,"richness,shannon_index,simpson_index\n");
  }
  // print header for averages
  fprintf(average_fp, "ProcNum,RcLine,Rep,Gen,MeanD,VarD,MeanH,VarH,MeanP,VarP,MeanM,VarM,N,1\n");
	
	for(rep=0; rep<ics.repeats; rep++) {
    marray_uint_set_zero(counts);
    marray_uint_set_zero(marginalized_driver_counts_per_generation);
    gsl_matrix_set_zero(dataM);
		// set initial N
		N = ics.Ninitial;
		// set all initial cells as unmutated
    indices[0] = 0;
    indices[1] = 0;
    indices[2] = 0;
    indices[3] = 0;
		marray_uint_set(counts, indices, N);
    min_occupied_D = 0;
    min_occupied_H = 0;
    min_occupied_P = 0;
    min_occupied_M = 0;
    max_occupied_D = 0;
    max_occupied_H = 0;
    max_occupied_P = 0;
    max_occupied_M = 0;
    time_first_sig_mut = ics.gens + 1;
    time_first_sig_ca = ics.gens + 1;

		for(i=0; i<ics.gens; i++) {

      N = do_cycle(gslr, counts, F, N, Md, Mh, Mp, Mm, 
                    &min_occupied_D, &min_occupied_H,
                    &min_occupied_P, &min_occupied_M,
                    &max_occupied_D, &max_occupied_H,
                    &max_occupied_P, &max_occupied_M,
                    &num_potentially_occupied_mutant_classes,
                    &potentially_occupied_mutant_classes,
                    index_of_potentially_occupied_mutant_class,
                    &mutant_class_probabilities,
                    &new_mutant_class_counts,
                    delta_counts, prob_of_cells_acquiring_num_muts_lo,
                    prob_of_cells_acquiring_num_muts_hi,
                    num_cells_acquiring_num_muts,
                    prob_that_mutation_falls_in_class, 
                    small_reciprocal_factorial,
                    reciprocal_num_remaining_loci,
                    indices, delta_indices,
                    new_indices);
      marginalized_indices[0] = i;
      if (ics.output_marginalized_counts==1) {
        compute_marginalized_driver_counts(counts, Md, 
                                  min_occupied_D, min_occupied_H,
                                  min_occupied_P, min_occupied_M,
                                  max_occupied_D, max_occupied_H, 
                                  max_occupied_P, max_occupied_M, indices,
                                  marginalized_driver_counts_per_generation, 
                                  marginalized_indices,
                                  marginalized_driver_count_fp);
        compute_marginalized_housekeeper_counts(counts, Mh, 
                                  min_occupied_D, min_occupied_H,
                                  min_occupied_P, min_occupied_M,
                                  max_occupied_D, max_occupied_H, 
                                  max_occupied_P, max_occupied_M, indices,
                                  marginalized_housekeeper_counts_per_generation, 
                                  marginalized_indices,
                                  marginalized_housekeeper_count_fp);
        compute_marginalized_passenger_counts(counts, Mp, 
                                  min_occupied_D, min_occupied_H,
                                  min_occupied_P, min_occupied_M,
                                  max_occupied_D, max_occupied_H, 
                                  max_occupied_P, max_occupied_M, indices,
                                  marginalized_passenger_counts_per_generation, 
                                  marginalized_indices,
                                  marginalized_passenger_count_fp);
        compute_marginalized_mutator_counts(counts, Mm, 
                                  min_occupied_D, min_occupied_H,
                                  min_occupied_P, min_occupied_M,
                                  max_occupied_D, max_occupied_H, 
                                  max_occupied_P, max_occupied_M, indices,
                                  marginalized_mutator_counts_per_generation, 
                                  marginalized_indices,
                                  marginalized_mutator_count_fp);
      }
			if(i%100==0) {
        report_data(process_number, counts, N, rep, i, 
                    min_occupied_D, min_occupied_H, min_occupied_P, min_occupied_M,
                    max_occupied_D, max_occupied_H, max_occupied_P, max_occupied_M,
                    time_first_sig_mut,time_first_sig_ca,indices,
                    &diversity_indices,
        fp);
				compute_averages(counts, N, 
                        min_occupied_D, min_occupied_H,
                        min_occupied_P, min_occupied_M,
                        max_occupied_D, max_occupied_H,
                        max_occupied_P, max_occupied_M, dataM, i/100, indices);
			}
			if(N>1e9) {
				printf("#process %d rep %d max tumor size reached\n", process_number,
              rep);
				break;
			}
      if(N==0) {
        printf("#process %d rep %d tumor went extinct\n", process_number, rep);
        break;
      }
			
			mutP = compute_number_of_mutator_cells(counts, 
                min_occupied_D, min_occupied_H, min_occupied_P, min_occupied_M, 
                max_occupied_D, max_occupied_H, max_occupied_P, max_occupied_M, 
                indices);
			prop = (double) mutP;
			prop /= (double) N;
			if(prop >= ics.mutProp) {
				printf("#process %d rep %d more than %g percent of tumor has a mutator phenotype\n",process_number, rep, prop*100);
				if(i < time_first_sig_mut) { 
					time_first_sig_mut = i;
				}
			}
			
			cancerP = compute_number_of_cells_classed_as_cancerous(counts, 
                                            min_occupied_D, min_occupied_H, 
                                            min_occupied_P, min_occupied_M, 
                                            max_occupied_D, max_occupied_H, 
                                            max_occupied_P, max_occupied_M, 
                                            indices);
			prop = (double) cancerP;
			prop /= (double) N;
			if(cancerP > 0) {
				if(prop >= ics.cancerProp) {
          if(i < time_first_sig_ca) { 
             time_first_sig_ca = i;
          }
					printf("#process %d rep %d made a cancer with at least %d driver mutations filling %g percent of tumor\n", process_number, rep, ics.mostagg,prop*100);
					break;
				}
			}
		}
    // print averages
    for(j=0; j<=GSL_MIN_INT(i,ics.gens-1)/100; j++) {
      denom = gsl_matrix_get(dataM, j, 9);
      reciprocal_denom = 1.0 / denom;
      fprintf(average_fp, "%d,%d,%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
              process_number,
              ics.run_number,
              rep,
              j*100,
              gsl_matrix_get(dataM, j, 0) * reciprocal_denom,
              gsl_matrix_get(dataM, j, 1) * reciprocal_denom,
              gsl_matrix_get(dataM, j, 2) * reciprocal_denom,
              gsl_matrix_get(dataM, j, 3) * reciprocal_denom,
              gsl_matrix_get(dataM, j, 4) * reciprocal_denom,
              gsl_matrix_get(dataM, j, 5) * reciprocal_denom,
              gsl_matrix_get(dataM, j, 6) * reciprocal_denom,
              gsl_matrix_get(dataM, j, 7) * reciprocal_denom,
              gsl_matrix_get(dataM, j, 8) * reciprocal_denom,
              denom);
    }
    if (time_first_sig_ca > ics.gens) {
      // Did not get to cancer
      fprintf(got_cancer_fp, "0,");
    } else {
      fprintf(got_cancer_fp, "1,");
    }
		report_data(process_number, counts, N, rep, i, 
                min_occupied_D, min_occupied_H, min_occupied_P, min_occupied_M,
                max_occupied_D, max_occupied_H, max_occupied_P, max_occupied_M,
                time_first_sig_mut,time_first_sig_ca,indices,
                &diversity_indices,
    got_cancer_fp);
	}
  fclose(got_cancer_fp);
  fclose(average_fp);
  if (ics.output_marginalized_counts==1) {
    fclose(marginalized_mutator_count_fp);
    fclose(marginalized_passenger_count_fp);
    fclose(marginalized_housekeeper_count_fp);
    fclose(marginalized_driver_count_fp);
  }
	fclose(fp);
	
	 
	// free memory
  free(new_mutant_class_counts);
  free(mutant_class_probabilities);
  free(potentially_occupied_mutant_classes);
  gsl_matrix_free(dataM);
	gsl_vector_uint_free(mutn);
	gsl_vector_uint_free(n);
  gsl_matrix_free(F);
  marray_uint_free(delta_counts);
  marray_uint_free(index_of_potentially_occupied_mutant_class);
  marray_uint_free(counts);
  gsl_rng_free(gslr);
	
  MPI_Finalize();
  return(0);
}

void compute_fitness_coeff(gsl_matrix *F, int Md, int Mh) {
	int i,j;
  double prev_entry, cur_entry, prev_driver, cur_driver;
  prev_driver = 1.0;
  prev_entry = prev_driver;
  gsl_matrix_set(F, 0, 0, 1.0);
	for(j=1; j<Mh+1; j++) {
    // Mutations to housekeeper genes decrease fitness
    cur_entry = prev_entry * (1-ics.sH);
    gsl_matrix_set(F, 0, j, cur_entry);
    prev_entry = cur_entry;
	}
	
  // Mutations to driver genes increase fitness
	for(i=1; i<Md+1; i++) {
    cur_driver = prev_driver * (1+ics.sD);
    gsl_matrix_set(F, i, 0, cur_driver);
    prev_driver = cur_driver;
    prev_entry = prev_driver;
    // Mutations to housekeeper genes decrease fitness
		for(j=1; j<Mh+1; j++) {
      cur_entry = prev_entry * (1-ics.sH);
      gsl_matrix_set(F, i, j, cur_entry);
      prev_entry = cur_entry;
		}
	}
}

void fill_powers(double *powers, double alpha, 
                          int highest_power) {
  double *cur_ptr;
  int i;
  double current_value;
  for (i = 0, cur_ptr = powers, current_value = 1.0; 
      i <= highest_power; 
      ++i, current_value *= alpha) {
    *cur_ptr++ = current_value; 
  }
}

void fill_prob_of_cells_acquiring_num_muts(
                                  double **prob_of_cells_acquiring_num_muts,
                                  int total_num_loci,
                                  double mu,
                                  double *mupower,
                                  double *oneminusmupower,
                                  double *small_reciprocal_factorial) {
  // prob_of_cells_acquiring_num_muts[i][j] is the probability of a cell
  // acquiring j new mutations, given that total_num_loci - i of the loci in
  // the cell have already been mutated
  // thus, there are i loci left to possibly mutate
  // this probability is 
  // (i choose j) mu ^ j * (1-mu) ^ (i - j)
  int i, j;
  double **prob_of_cells_acquiring_zero_muts_ptr;
  double *prob_of_cells_acquiring_num_muts_ptr;
  double *mu_power_ptr;
  double *top_oneminusmu_power_ptr;
  double *oneminusmu_power_ptr;
  double *small_reciprocal_factorial_ptr;
  int i_choose_j_numerator;
  fill_powers(mupower, mu, MAXMUTSPERCELL);
  fill_powers(oneminusmupower, 1.0 - mu, total_num_loci);
  for (i = 0, 
      prob_of_cells_acquiring_zero_muts_ptr = prob_of_cells_acquiring_num_muts,
      top_oneminusmu_power_ptr = oneminusmupower;
      i <= total_num_loci; 
      ++i,
      ++prob_of_cells_acquiring_zero_muts_ptr,
      ++top_oneminusmu_power_ptr) {
    // (i choose 0) is just 1
    i_choose_j_numerator = 1; 
    for (j = 0,
        prob_of_cells_acquiring_num_muts_ptr =
        *prob_of_cells_acquiring_zero_muts_ptr,
        mu_power_ptr = mupower,
        oneminusmu_power_ptr = top_oneminusmu_power_ptr,
        small_reciprocal_factorial_ptr = small_reciprocal_factorial;
        j <= MAXMUTSPERCELL;
        ++j,
        ++prob_of_cells_acquiring_num_muts_ptr,
        ++mu_power_ptr,
        --oneminusmu_power_ptr, 
        ++small_reciprocal_factorial_ptr) {
      if (j > i) {
        *prob_of_cells_acquiring_num_muts_ptr = 0.0;
      } else {
        *prob_of_cells_acquiring_num_muts_ptr = 
          *mu_power_ptr *
          *oneminusmu_power_ptr *
          i_choose_j_numerator *
          *small_reciprocal_factorial_ptr;
      }
      // This is for use in the next iteration of this loop
      i_choose_j_numerator *= (i - j);
    }
  }
}

void fill_reciprocal_lookup(double *reciprocal_num_remaining_loci, 
                            int total_num_loci) {
  int i;
  double *cur_ptr;
  for (i = 0, cur_ptr = reciprocal_num_remaining_loci;
      i <= total_num_loci;
      ++i) {
    *cur_ptr++ = 1.0 / i;
  }
}

double compute_mean_fitness(marray_uint *counts, gsl_matrix *F, 
                      int min_occupied_D, int min_occupied_H,
                      int min_occupied_P, int min_occupied_M,
                      int max_occupied_D, int max_occupied_H,
                      int max_occupied_P, int max_occupied_M,
                      double reciprocal_population_size, size_t *indices) {
  double fitness, mean_fitness;
  mean_fitness = 0.0;
  for (indices[0] = min_occupied_D; indices[0] <= max_occupied_D; indices[0]++) {
    for (indices[1] = min_occupied_H; indices[1] <= max_occupied_H; indices[1]++) {
      // This fitness applies to each mutant class with the specific number of
      // driver mutations indices[0] and the specific number of housekeeper
      // mutations indices[1].
      fitness = gsl_matrix_get(F, indices[0], indices[1]);
      for (indices[2] = min_occupied_P; indices[2] <= max_occupied_P; indices[2]++) {
        for (indices[3] = min_occupied_M; indices[3] <= max_occupied_M; indices[3]++) {
          mean_fitness += fitness * marray_uint_get(counts, indices);
        }
      }
    }
  }
  mean_fitness *= reciprocal_population_size;
  return mean_fitness;
}

double compute_w_bar(marray_uint *counts, gsl_matrix *F, 
                      int min_occupied_D, int min_occupied_H,
                      int min_occupied_P, int min_occupied_M,
                      int max_occupied_D, int max_occupied_H,
                      int max_occupied_P, int max_occupied_M,
                      double reciprocal_population_size,
                      double reciprocal_mean_fitness, size_t *indices) {
  double fitness, w_bar;
  w_bar = 0.0;
  for (indices[0] = min_occupied_D; indices[0] <= max_occupied_D; indices[0]++) {
    for (indices[1] = min_occupied_H; indices[1] <= max_occupied_H; indices[1]++) {
      // This fitness applies to each mutant class with the specific number of
      // driver mutations indices[0] and the specific number of housekeeper
      // mutations indices[1].
      fitness = gsl_matrix_get(F, indices[0], indices[1]);
      for (indices[2] = min_occupied_P; indices[2] <= max_occupied_P; indices[2]++) {
        for (indices[3] = min_occupied_M; indices[3] <= max_occupied_M; indices[3]++) {
          w_bar += fitness * marray_uint_get(counts, indices);
        }
      }
    }
  }
  w_bar *= reciprocal_mean_fitness;
  w_bar *= reciprocal_population_size;
  return w_bar;
}

unsigned int do_cycle(gsl_rng *gslr, marray_uint *counts, gsl_matrix *F,
    unsigned int N, int Md, int Mh, int Mp, int Mm, 
    int *min_occupied_D_ptr, int *min_occupied_H_ptr,
    int *min_occupied_P_ptr, int *min_occupied_M_ptr,
    int *max_occupied_D_ptr, int *max_occupied_H_ptr,
    int *max_occupied_P_ptr, int *max_occupied_M_ptr,
    int *num_potentially_occupied_mutant_classes_ptr,
    mutant_class_indices **potentially_occupied_mutant_classes_ptr,
    marray_uint *index_of_potentially_occupied_mutant_class,
    double **mutant_class_probabilities_ptr,
    unsigned int **new_mutant_class_counts_ptr,
    marray_uint *delta_counts,
    double *prob_of_cells_acquiring_num_muts_lo[],
    double *prob_of_cells_acquiring_num_muts_hi[],
    unsigned int num_cells_acquiring_num_muts[],
    double prob_that_mutation_falls_in_class[],
    double *small_reciprocal_factorial,
    double *reciprocal_num_remaining_loci,
    size_t *indices,
    size_t *delta_indices, 
    size_t *new_indices) {
  double mean_fitness;
  double reciprocal_population_size;
	
  switch (ics.growthModel) {
    case 'I':
      // Here, there is no constraint on population size, and each mutant class
      // independently (I) grows exponentially
      // update total counts in each class using corresponding fitness
      N = grow_total_counts(gslr, counts, F, *max_occupied_D_ptr,
                        *max_occupied_H_ptr, *max_occupied_P_ptr,
                        *max_occupied_M_ptr, indices);
      // Mutate the cells
      mutate_counts(gslr, counts, Md, Mh, Mp, Mm, 
        max_occupied_D_ptr, max_occupied_H_ptr,
        max_occupied_P_ptr, max_occupied_M_ptr,
        delta_counts,
        prob_of_cells_acquiring_num_muts_lo,
        prob_of_cells_acquiring_num_muts_hi,
        num_cells_acquiring_num_muts,
        prob_that_mutation_falls_in_class, 
        indices, delta_indices, new_indices);
      return(N);
    case 'E':
      // Here, the population size is constrained to grow exponentially.  This
      // model is here primarily for comparison with other papers.
      // N = exp(Ct) => dN/dt = C exp(Ct) = C N
      N = gsl_ran_poisson(gslr, N + (unsigned int)rint(ics.beta * N));
      break;
    case 'B':
      // Here, the population size grows exponentially according to the
      // absolute fitness of the total population, and this population size is
      // also taken as the carrying capacity constraint, within which the
      // clones compete according to their relative fitnesses (the same ones
      // that were used to compute the absolute fitness).  This model is here
      // to duplicate the one in the original Beerenwinkel paper.
      reciprocal_population_size = 1.0 / N;
      mean_fitness = compute_mean_fitness(counts, F,
                                                          *min_occupied_D_ptr,
                                                          *max_occupied_D_ptr,
                                                          *min_occupied_H_ptr,
                                                          *max_occupied_H_ptr,
                                                          *min_occupied_P_ptr,
                                                          *max_occupied_P_ptr,
                                                          *min_occupied_M_ptr,
                                                          *max_occupied_M_ptr,
                                                          reciprocal_population_size,
                                                          indices);
      N = gsl_ran_poisson(gslr, N + (unsigned int)rint(ics.beta * mean_fitness * N));
    case 'C':
      // cubic population growth
      // N = Ct^3 => dN/dt = 3Ct^2 = K N^(2/3)
      // Think of cells expanding volumetrically into a spatial territory: the
      // cells at the periphery (the surface of the ball of cells) are the ones
      // able to expand the available territory.
      N = gsl_ran_poisson(gslr, 
                          N + (unsigned int)rint(ics.beta * pow(cbrt(N),2)));
      break;
    case 'Q':
      // quadratic population growth
      // N = Ct^2 => dN/dt = 2Ct = K N^(1/2)
      // Think of cells constrained to expand along a surface (e.g., by the
      // basement membrane): the cells at the periphery (i.e., along the
      // perimeter) are the ones able to expand the available territory.
      N = gsl_ran_poisson(gslr, 
                          N + (unsigned int)rint(ics.beta * sqrt(N)));
      break;
    case 'K':
      // constant population size
      // N = C
      // The total number of cells is constrained by the availability of some
      // resource which does not grow.
      N = gsl_ran_poisson(gslr, N);
    default:
      break;
  }
  update_class_counts(gslr, N, counts, F, Md, Mh, Mp, Mm,
                      min_occupied_D_ptr, min_occupied_H_ptr,
                      min_occupied_P_ptr, min_occupied_M_ptr,
                      max_occupied_D_ptr, max_occupied_H_ptr,
                      max_occupied_P_ptr, max_occupied_M_ptr,
                      num_potentially_occupied_mutant_classes_ptr,
                      potentially_occupied_mutant_classes_ptr,
                      index_of_potentially_occupied_mutant_class,
                      mutant_class_probabilities_ptr,
                      new_mutant_class_counts_ptr,
                      prob_of_cells_acquiring_num_muts_lo,
                      prob_of_cells_acquiring_num_muts_hi,
                      small_reciprocal_factorial,
                      reciprocal_num_remaining_loci,
                      indices,
                      delta_indices, new_indices);
	
	return(N);
}

unsigned int grow_total_counts(gsl_rng *gslr, marray_uint *counts, gsl_matrix *F, int Md, int Mh, int Mp, int Mm, size_t *indices) {
  double fitness;
  unsigned int count, new_count;
  unsigned int newN;
  newN = 0;
  for (indices[0] = 0; indices[0] <= Md; indices[0]++) {
    for (indices[1] = 0; indices[1] <= Mh; indices[1]++) {
      fitness = gsl_matrix_get(F, indices[0], indices[1]);
      for (indices[2] = 0; indices[2] <= Mp; indices[2]++) {
        for (indices[3] = 0; indices[3] <= Mm; indices[3]++) {
          count = marray_uint_get(counts, indices);
          if (count > 0) {
            new_count = gsl_ran_poisson(gslr, rint(count * fitness * ics.beta));
            newN += new_count;
            marray_uint_set(counts, indices, new_count);
          }
        }
      }
    }
  }
  return(newN);
}

void mutate_counts(gsl_rng *gslr, marray_uint *counts, int Md, int Mh, int Mp, int Mm, 
    int *max_occupied_D_ptr, int *max_occupied_H_ptr, int *max_occupied_P_ptr,
    int *max_occupied_M_ptr,
    marray_uint *delta_counts,
    double *prob_of_cells_acquiring_num_muts_lo[],
    double *prob_of_cells_acquiring_num_muts_hi[],
    unsigned int num_cells_acquiring_num_muts[],
    double prob_that_mutation_falls_in_class[], 
    size_t *indices, size_t *delta_indices, size_t *new_indices) {
  // We traverse the mutation classes in an order analogous to Cantor's
  // diagonalization process, so that when we update the mutation class 
  // (i, j, k, l), we have already updated all mutation classes 
  // (i', j', k', l') with a greater total number of mutations.  This is
  // because when we introduce new mutations in a cell of class (i, j, k, l),
  // we shift this cell from class (i, j, k, l) to some other class 
  // (i', j', k', l') with a greater total number of mutations, incrementing
  // the count in that class.  If we hadn't already updated mutation class 
  // (i', j', k', l'), then this cell could get mutated again when we came to
  // update class (i', j', k', l'), which would not be a correct implementation
  // of our model.  
  int max_total_mutations;
  int maxD, maxH, maxP, maxM;
  int total_mutations, total_nondriver_mutations, total_nondriver_nonhousekeeper_mutations;
  int max_num_mutator_mutations;
  maxD = GSL_MAX_INT(*max_occupied_D_ptr + MAXMUTSPERCELL, Md);
  maxH = GSL_MAX_INT(*max_occupied_H_ptr + MAXMUTSPERCELL, Mh);
  maxP = GSL_MAX_INT(*max_occupied_P_ptr + MAXMUTSPERCELL, Mp);
  maxM = GSL_MAX_INT(*max_occupied_M_ptr + MAXMUTSPERCELL, Mm);
  max_total_mutations = maxD + maxH + maxP + maxM;
  for (total_mutations = max_total_mutations; total_mutations >= 0; total_mutations--) {
    for (indices[0] = GSL_MIN_INT(total_mutations, maxD), 
        total_nondriver_mutations = total_mutations - indices[0];         // indices[0] + total_nondriver_mutations = total_mutations
        indices[0] >= 0 && total_nondriver_mutations <= total_mutations;
        --indices[0], ++total_nondriver_mutations) {                      // the above remains invariant
      for (indices[1] = GSL_MIN_INT(total_nondriver_mutations, maxH),
          total_nondriver_nonhousekeeper_mutations = total_nondriver_mutations - indices[1];  // indices[1] + total_nondriver_nonhousekeeper_mutations = total_nondriver_mutations
          indices[1] >= 0 && total_nondriver_nonhousekeeper_mutations <= total_nondriver_mutations;
          --indices[1], ++total_nondriver_nonhousekeeper_mutations) {                         // the above remains invariant
        max_num_mutator_mutations = GSL_MIN_INT(total_nondriver_nonhousekeeper_mutations, 
                                                maxM);
        for (indices[2] = GSL_MIN_INT(total_nondriver_nonhousekeeper_mutations, 
              maxP),
            indices[3] = total_nondriver_nonhousekeeper_mutations - indices[2];                     // indices[2] + indices[3] = total_nondriver_nonhousekeeper_mutations
            indices[2] >= 0 && indices[3] <= max_num_mutator_mutations;
            --indices[2], ++indices[3]) {                                                           // the above remains invariant
          // Substituting in all the invariants, we have:
          // indices[0] + indices[1] + indices[2] + indices[3] = total_mutations
          // as well as
          // 0 <= indices[0] <= maxD
          // 0 <= indices[1] <= maxH
          // 0 <= indices[2] <= maxP
          // 0 <= indices[3] <= maxM
          // at every iteration of the loop.
          if (indices[3] >= ics.mut_max) {
            // The cells in this class have enough mutator mutations to have
            // the mutator phenotype
            mutate_counts_in_class(gslr, counts, Md, Mh, Mp, Mm, 
                                  max_occupied_D_ptr, max_occupied_H_ptr,
                                  max_occupied_P_ptr, max_occupied_M_ptr,
                                  delta_counts, 
                                  prob_of_cells_acquiring_num_muts_hi,
                                  num_cells_acquiring_num_muts,
                                  prob_that_mutation_falls_in_class,
                                  indices, delta_indices, new_indices, 1);
          } else {
            // The cells in this class do not have the mutator phenotype
            mutate_counts_in_class(gslr, counts, Md, Mh, Mp, Mm, 
                                  max_occupied_D_ptr, max_occupied_H_ptr,
                                  max_occupied_P_ptr, max_occupied_M_ptr,
                                  delta_counts, 
                                  prob_of_cells_acquiring_num_muts_lo,
                                  num_cells_acquiring_num_muts,
                                  prob_that_mutation_falls_in_class,
                                  indices, delta_indices, new_indices, 0);
          }
        }
      }
    }
  }
}

void mutate_counts_in_class(gsl_rng *gslr, marray_uint *counts, 
    int Md, int Mh, int Mp, int Mm, 
    int *max_occupied_D_ptr, int *max_occupied_H_ptr,
    int *max_occupied_P_ptr, int *max_occupied_M_ptr,
    marray_uint *delta_counts,
    double *prob_of_cells_acquiring_num_muts[],
    unsigned int num_cells_acquiring_num_muts[],
    double prob_that_mutation_falls_in_class[], 
    size_t *indices, size_t *delta_indices, size_t *new_indices, 
    int mutation_rate_level) {
  int num_mutations;
  unsigned int Nclass;
  unsigned int delta_count;
  bool couldBeOccupied;
  Nclass = marray_uint_get(counts, indices);
  if (Nclass == 0) {
    return;
  }
  // We remove these Nclass cells; we will redistribute them, after mutating
  marray_uint_set(counts, indices, 0);
  // Sample the number of cells acquiring each possible number of mutations
  // from the appropriate multinomial distribution
  gsl_ran_multinomial(gslr, MAXMUTSPERCELL+mutation_rate_level, Nclass,
          *prob_of_cells_acquiring_num_muts, num_cells_acquiring_num_muts);
  marray_uint_set_zero(delta_counts);
  for (num_mutations = MAXMUTSPERCELL-1+mutation_rate_level; num_mutations >= 0; num_mutations--) {
    // Make the appropriate number of cells acquire num_mutations
    delta_indices[0] = 0;
    delta_indices[1] = 0;
    delta_indices[2] = 0;
    delta_indices[3] = 0;
    distribute_mutations(gslr, delta_counts, delta_indices,
      prob_that_mutation_falls_in_class,
      num_cells_acquiring_num_muts[num_mutations], num_mutations);
  }
  // Update counts by moving cells into new mutation classes
  // If a cell would have more than the maximum number of mutations of a given
  // class, give it the maximum number of mutations for that class
  for (delta_indices[0] = 0; 
      delta_indices[0] < MAXMUTSPERCELL+mutation_rate_level;
      delta_indices[0]++) {
    new_indices[0] = GSL_MIN_INT(indices[0] + delta_indices[0], Md);
    for (delta_indices[1] = 0; 
        delta_indices[1] < MAXMUTSPERCELL+mutation_rate_level;
        delta_indices[1]++) {
      new_indices[1] = GSL_MIN_INT(indices[1] + delta_indices[1], Mh);
      for (delta_indices[2] = 0; 
          delta_indices[2] < MAXMUTSPERCELL+mutation_rate_level;
          delta_indices[2]++) {
        new_indices[2] = GSL_MIN_INT(indices[2] + delta_indices[2], Mp);
        for (delta_indices[3] = 0; 
            delta_indices[3] < MAXMUTSPERCELL+mutation_rate_level;
            delta_indices[3]++) {
          delta_count = marray_uint_get(delta_counts, delta_indices);
          if (delta_count > 0) {
            new_indices[3] = GSL_MIN_INT(indices[3] + delta_indices[3], Mm);
            couldBeOccupied = true;
            if (new_indices[0] > *max_occupied_D_ptr) {
              *max_occupied_D_ptr = new_indices[0];
              couldBeOccupied = false;
            }
            if (new_indices[1] > *max_occupied_H_ptr) {
              *max_occupied_H_ptr = new_indices[1];
              couldBeOccupied = false;
            }
            if (new_indices[2] > *max_occupied_P_ptr) {
              *max_occupied_P_ptr = new_indices[2];
              couldBeOccupied = false;
            }
            if (new_indices[3] > *max_occupied_M_ptr) {
              *max_occupied_M_ptr = new_indices[3]; 
              couldBeOccupied = false;
            }
            if (couldBeOccupied) {
              marray_uint_set(counts, new_indices, 
                              marray_uint_get(counts, new_indices) + delta_count);
            } else {
              marray_uint_set(counts, new_indices, delta_count);
            }
          }
        }
      }
    }
  }
}

void distribute_mutations(gsl_rng *gslr, marray_uint *delta_counts, size_t *delta_indices,
                          double *prob_that_mutation_falls_in_class, unsigned int num_cells, 
                          unsigned int num_mutations_left) {
  if (num_cells == 0) {
    return;
  }
  if (num_mutations_left == 0) {
    marray_uint_set(delta_counts, delta_indices, marray_uint_get(delta_counts, delta_indices) + num_cells);
    return;
  }
  unsigned int num_cells_acquiring_muts_in_class[4];
  int class;
  // Sample the number of cells acquiring a new mutation to a class of genes
  // from the multinomial distribution
  gsl_ran_multinomial(gslr, 4, num_cells, prob_that_mutation_falls_in_class, num_cells_acquiring_muts_in_class);

  for (class = 0; class < 4; ++class) {
    delta_indices[class] += 1;
    distribute_mutations(gslr, delta_counts, delta_indices, prob_that_mutation_falls_in_class,
                          num_cells_acquiring_muts_in_class[class], num_mutations_left-1);
    delta_indices[class] -= 1;
  }
  return;
}

void update_class_counts(gsl_rng *gslr, int N, marray_uint *counts, 
                        gsl_matrix *F,
                        int Md, int Mh, int Mp, int Mm,
                        int *min_occupied_D_ptr, int *min_occupied_H_ptr,
                        int *min_occupied_P_ptr, int *min_occupied_M_ptr,
                        int *max_occupied_D_ptr, int *max_occupied_H_ptr,
                        int *max_occupied_P_ptr, int *max_occupied_M_ptr,
                        int *num_potentially_occupied_mutant_classes_ptr,
                        // These are the 4-dimensional indices of each of the
                        // occupied classes.  The 4-dimensional indices index
                        // into the 4-dimension counts marray.  This array
                        // itself is linear.
                        mutant_class_indices **potentially_occupied_mutant_classes_ptr,
                        // These are the linear indices of each of the occupied
                        // classes.  The linear indices index into the
                        // new_mutant_class_counts_ptr array.  This marray is
                        // indexed by 4-dimensional indices.
                        marray_uint *index_of_potentially_occupied_mutant_class,
                        // These are the probabilities of each of the occupied
                        // classes in the next generation.  This array is
                        // linearly indexed.
                        double **mutant_class_probabilities_ptr,
                        // These are the counts of each class in the next
                        // generation, sampled from the multinomial probability
                        // distribution.  This array is linearly indexed.
                        unsigned int **new_mutant_class_counts_ptr,
                        double *prob_of_cells_acquiring_num_muts_lo[],
                        double *prob_of_cells_acquiring_num_muts_hi[],
                        double *small_reciprocal_factorial,
                        double *reciprocal_num_remaining_loci,
                        size_t *indices, size_t *delta_indices, 
                        size_t *new_indices) {
  int maxD, maxH, maxP, maxM, num_assigned_indices, class, total_num_loci;
  int maximum_of_class[4];
  unsigned int idx;
  int num_potentially_occupied_mutant_classes;
  size_t *index_ptr, *index_ptr2, *index_ptr3;
  mutant_class_indices *class_indices_ptr;
  double *class_probability_ptr;
  unsigned int *new_class_count_ptr;
  unsigned int count;
  double fitness;
  int mutation_number, num_new_mutations;
  int mutation_level;
  unsigned int mutation_code, max_mutation_code;
  unsigned int decoded_mutation_code;
  double base_probability_increment, new_class_probability_increment;
  double mutation_combination_factor;
  double Z; // Normalization constant for the mutant class probabilities
  double reciprocal_Z;
  mutant_class_indices *potentially_occupied_mutant_classes;
  double *mutant_class_probabilities;
  unsigned int *new_mutant_class_counts;
  int num_mutated_loci;
  double *prob_of_cells_acquiring_num_muts_ptr;
  double *reciprocal_cur_num_remaining_loci_ptr;
  maximum_of_class[0] = Md;
  maximum_of_class[1] = Mh;
  maximum_of_class[2] = Mp;
  maximum_of_class[3] = Mm;
  total_num_loci = (Md + 1) + (Mh + 1) + (Mp + 1) + (Mm + 1);
  maxD = GSL_MIN_INT(*max_occupied_D_ptr + MAXMUTSPERCELL, Md);
  maxH = GSL_MIN_INT(*max_occupied_H_ptr + MAXMUTSPERCELL, Mh);
  maxP = GSL_MIN_INT(*max_occupied_P_ptr + MAXMUTSPERCELL, Mp);
  maxM = GSL_MIN_INT(*max_occupied_M_ptr + MAXMUTSPERCELL, Mm);
  // We overestimate the number of potentially occupied mutant classes by
  // making a rectilinear box in the grid covering the current occupied mutant
  // classes, and then adding MAXMUTSPERCELL in each direction.
  // Because we expect travelling waves, we expect this number will be much
  // smaller than the total number of mutant classes.
  num_potentially_occupied_mutant_classes 
    = (maxD - *min_occupied_D_ptr + 1) 
    * (maxH - *min_occupied_H_ptr + 1) 
    * (maxP - *min_occupied_P_ptr + 1) 
    * (maxM - *min_occupied_M_ptr + 1);
  if (num_potentially_occupied_mutant_classes >
          *num_potentially_occupied_mutant_classes_ptr) {
    // The current arrays may not fit all the potentially occupied mutant
    // classes.  Keep adding CLASSBLOCKSIZE to the number of entries available
    // until there are enough for all the potentially occupied mutant classes.
    while (num_potentially_occupied_mutant_classes >
            *num_potentially_occupied_mutant_classes_ptr) {
      *num_potentially_occupied_mutant_classes_ptr += CLASSBLOCKSIZE;
    }
    // Now extend each of the arrays to have at least
    // *num_potentially_occupied_mutant_classes_ptr entries.
    *potentially_occupied_mutant_classes_ptr 
      = realloc(*potentially_occupied_mutant_classes_ptr,
                *num_potentially_occupied_mutant_classes_ptr 
                  * sizeof(mutant_class_indices));
    if (!(*potentially_occupied_mutant_classes_ptr)) {
      printf("Unable to reallocate memory for potentially_occupied_mutant_classes\n");
      exit(1);
    }
    *mutant_class_probabilities_ptr = realloc(*mutant_class_probabilities_ptr,
        *num_potentially_occupied_mutant_classes_ptr * sizeof(double));
    if (!(*mutant_class_probabilities_ptr)) {
      printf("Unable to reallocate memory for mutant_class_probabilities\n");
      exit(1);
    }
    *new_mutant_class_counts_ptr = realloc(*new_mutant_class_counts_ptr,
        *num_potentially_occupied_mutant_classes_ptr * sizeof(unsigned int));
    if (!(*new_mutant_class_counts_ptr)) {
      printf("Unable to reallocate memory for new_mutant_class_counts\n");
      exit(1);
    }
  }
  potentially_occupied_mutant_classes 
    = *potentially_occupied_mutant_classes_ptr;
  mutant_class_probabilities = *mutant_class_probabilities_ptr;
  new_mutant_class_counts = *new_mutant_class_counts_ptr;
  // idx is a linear index into each of these arrays
  // initialize each of the linear arrays to zero
  for (idx = 0, class_indices_ptr = potentially_occupied_mutant_classes,
      class_probability_ptr = mutant_class_probabilities,
      new_class_count_ptr = new_mutant_class_counts;
      idx < *num_potentially_occupied_mutant_classes_ptr;
      ++idx,
      ++class_indices_ptr,
      ++class_probability_ptr,
      ++new_class_count_ptr) {
    index_ptr = class_indices_ptr->indices;
    *index_ptr++ = 0;
    *index_ptr++ = 0;
    *index_ptr++ = 0;
    *index_ptr = 0;
    *class_probability_ptr = 0.0;
    *new_class_count_ptr = 0;
  }
  num_assigned_indices = 0;
  // class_indices_ptr will keep moving through the linear array, always
  // pointing to the next unassigned entry
  class_indices_ptr = potentially_occupied_mutant_classes;
  Z = 0.0;
  // indices[0], indices[1], indices[2], and indices[3] correspond to p, q, r,
  // u in the formula in the paper;
  // new_indices[0], new_indices[1], new_indices[2], new_indices[3] correspond
  // to i, j, k, l in the formula in the paper
  // Here, the currently occupied mutant class has (p, q, r, u) driver,
  // housekeeper, passenger, and mutator mutations respectively.  The
  // destination mutant class has (i, j, k, l) respective mutations.
  for (indices[0] = *min_occupied_D_ptr; indices[0] <= *max_occupied_D_ptr; indices[0]++) {
    for (indices[1] = *min_occupied_H_ptr; indices[1] <= *max_occupied_H_ptr; indices[1]++) {
      // This fitness applies to each mutant class with the specific number of
      // driver mutations indices[0] and the specific number of housekeeper
      // mutations indices[1].
      fitness = gsl_matrix_get(F, indices[0], indices[1]);
      for (indices[2] = *min_occupied_P_ptr; indices[2] <= *max_occupied_P_ptr; indices[2]++) {
        for (indices[3] = *min_occupied_M_ptr; indices[3] <= *max_occupied_M_ptr; indices[3]++) {
          count = marray_uint_get(counts, indices);
          if (count > 0) {  // this class is occupied
            if (indices[3] >= ics.mut_max) {
              // ics.mut_max is the number of mutations needed to switch to a
              // mutator phenotype
              mutation_level = 1;
            } else {
              mutation_level = 0;
            }
            // First do the case of no mutations
            // set idx to the index in the linear array corresponding to this
            // mutant class
            // At the first iteration of the outer loop,
            // index_of_potentially_occupied_mutant_class is all zero, so the
            // returned idx will always be zero.
            // As we go through the outer loop, we come across potentially
            // occupied mutant classes.  If they have not already been mapped
            // to a position in the linear array, we add a new entry at the end
            // of the array and assign it the indices of this potentially
            // occupied mutant class.  We update the mapping between the
            // 4-dimensional indices and the index into the linear array.
            // The end of the linear array is given by num_assigned_indices.
            idx = marray_uint_get(index_of_potentially_occupied_mutant_class,
                                          indices);
            if (idx == 0) {
              // Possibly no index exists in the linear array corresponding to
              // this mutant class, so we need to add a new class to the linear
              // array.  (It also may be that the idx really is zero, but we'll
              // just have to do that case again.)
              marray_uint_set(index_of_potentially_occupied_mutant_class,
                              indices, num_assigned_indices);
              // class_indices_ptr points to the end of the linear array
              // potentially_occupied_mutant_classes
              // Assign indices to class_indices_ptr->indices
              index_ptr = class_indices_ptr->indices;
              index_ptr2 = indices;
              *index_ptr++ = *index_ptr2++;
              *index_ptr++ = *index_ptr2++;
              *index_ptr++ = *index_ptr2++;
              *index_ptr = *index_ptr2;
              idx = num_assigned_indices;
              ++class_indices_ptr;
              ++num_assigned_indices;
            }
            // One of the mutually exclusive ways of getting a cell in this
            // mutant class, is from a previous generation cell of this class
            // replicating and not mutating.
            // prob_of_cells_acquiring_num_muts[i][j] is the probability of a cell
            // acquiring j new mutations, given that total_num_loci - i of the loci in
            // the cell have already been mutated
            num_mutated_loci 
              = indices[0] + indices[1] + indices[2] + indices[3];
            reciprocal_cur_num_remaining_loci_ptr 
              = &(reciprocal_num_remaining_loci[total_num_loci -
                                                num_mutated_loci]);
            if (mutation_level == 0) {
              prob_of_cells_acquiring_num_muts_ptr =
                prob_of_cells_acquiring_num_muts_lo[total_num_loci -
                  num_mutated_loci];
            } else {
              prob_of_cells_acquiring_num_muts_ptr =
                prob_of_cells_acquiring_num_muts_hi[total_num_loci -
                  num_mutated_loci];
            }
            new_class_probability_increment =
              count       // proportion of mutant class in previous generation
              * fitness   // fitness of this mutant class
              * *prob_of_cells_acquiring_num_muts_ptr++  // prob of no mutations
              ;
            mutant_class_probabilities[idx] += new_class_probability_increment;
            Z += new_class_probability_increment;
            max_mutation_code = 1;
            // Loop through all possible ways of getting num_new_mutations
            // The combinations of mutations are encoded in a base-4 number
            // with num_new_mutations digits
            // Each digit signifies a new mutation of type:
            // 0 - driver
            // 1 - housekeeper
            // 2 - passenger
            // 3 - mutant
            // For instance, the base-4 number 203 indicates a new passenger
            // mutation, a new driver mutation, and a new mutator mutation.
            // By looping through the mutation codes, we will automatically
            // take into account the number of different ways of getting a
            // particular combination of new mutations of different classes.
            for (num_new_mutations = 1; 
                  num_new_mutations <= MAXMUTSPERCELL - 1 + mutation_level;
                  num_new_mutations++) {
                base_probability_increment =
                  count       // proportion of mutant class in previous generation
                  * fitness   // fitness of this mutant class
                  * *prob_of_cells_acquiring_num_muts_ptr++;
              max_mutation_code *= 4;
              for (mutation_code = 0; mutation_code < max_mutation_code;
                  ++mutation_code) {
                mutation_combination_factor = 1.0;
                // Zero out delta_indices
                // delta_indices tells us where is the destination mutant class
                // compared to the current (occupied) mutant class
                index_ptr = delta_indices;
                *index_ptr++ = 0;
                *index_ptr++ = 0;
                *index_ptr++ = 0;
                *index_ptr = 0;
                decoded_mutation_code = mutation_code;
                // Read off the base-4 digits of the mutation code using
                // decoded_mutation_code
                for (mutation_number = 0; mutation_number < num_new_mutations;
                      ++mutation_number) {
                  // The least significant base-4 digit tells us the type of
                  // the current mutation
                  class = decoded_mutation_code & 3;
                  // We multiply by the proportion of the remaining loci that
                  // are of this type:

                  // Multiply by the number of previously unmutated loci of
                  // this type in which this mutation could have occurred
                  mutation_combination_factor *= GSL_MAX_INT(0, (maximum_of_class[class] -
                                                  (indices[class]+delta_indices[class])));
                  if (mutation_combination_factor == 0.0)
                    break;
                  // Divide by the total number of previously unmutated loci
                  mutation_combination_factor 
                    *= *reciprocal_cur_num_remaining_loci_ptr;
                  delta_indices[class] += 1;
                  // Shift so the next least significant base-4 digit ends up in
                  // the least significant place
                  decoded_mutation_code >>= 2;
                }
                if (mutation_combination_factor != 0.0) {
                  for (class = 0; class < 4; ++class) {
                    // We don't care which order the loci in the class were
                    // mutated, so divide by (# of newly mutated loci in class)!
                    mutation_combination_factor 
                      *= small_reciprocal_factorial[delta_indices[class]];
                  }
                  // Assign indices + delta_indices to new_indices
                  index_ptr = new_indices;
                  index_ptr2 = indices;
                  index_ptr3 = delta_indices;
                  *index_ptr++ = *index_ptr2++ + *index_ptr3++;
                  *index_ptr++ = *index_ptr2++ + *index_ptr3++;
                  *index_ptr++ = *index_ptr2++ + *index_ptr3++;
                  *index_ptr = *index_ptr2 + *index_ptr3;
                  idx = marray_uint_get(index_of_potentially_occupied_mutant_class,
                                                new_indices);
                  if (idx == 0) {
                    marray_uint_set(index_of_potentially_occupied_mutant_class,
                                    new_indices, num_assigned_indices);
                    // Assign new_indices to class_indices_ptr
                    index_ptr = class_indices_ptr->indices;
                    index_ptr2 = new_indices;
                    *index_ptr++ = *index_ptr2++;
                    *index_ptr++ = *index_ptr2++;
                    *index_ptr++ = *index_ptr2++;
                    *index_ptr = *index_ptr2;
                    idx = num_assigned_indices;
                    ++class_indices_ptr;
                    ++num_assigned_indices;
                  }
                  new_class_probability_increment 
                    = base_probability_increment * mutation_combination_factor;
                  mutant_class_probabilities[idx] 
                    += new_class_probability_increment;
                  Z += new_class_probability_increment;
                }
              }
              // We have mutated one locus, so num_mutated_loci has increased
              // by one, and the number of remaining loci has decreased by one.
              --reciprocal_cur_num_remaining_loci_ptr;
            }
          }
        }
      }
    }
  }
  // Normalize the probabilities
  reciprocal_Z = 1.0 / Z;
  for (idx = 0, class_probability_ptr = mutant_class_probabilities; 
        idx < num_assigned_indices; ++idx, ++class_probability_ptr) {
    *class_probability_ptr *= reciprocal_Z;
  }
  // Sample the multinomial distribution
  gsl_ran_multinomial(gslr, num_assigned_indices, N,
          mutant_class_probabilities, new_mutant_class_counts);
  // Assign these back into counts
  *min_occupied_D_ptr = Md;
  *max_occupied_D_ptr = 0;
  *min_occupied_H_ptr = Mh;
  *max_occupied_H_ptr = 0;
  *min_occupied_P_ptr = Mp;
  *max_occupied_P_ptr = 0;
  *min_occupied_M_ptr = Mm;
  *max_occupied_M_ptr = 0;
  for (idx = 0, class_indices_ptr = potentially_occupied_mutant_classes,
      new_class_count_ptr = new_mutant_class_counts;
      idx < num_assigned_indices;
      ++idx,
      ++class_indices_ptr,
      ++new_class_count_ptr) {
    marray_uint_set(counts, class_indices_ptr->indices, *new_class_count_ptr);
    if (*new_class_count_ptr > 0) {
      index_ptr = class_indices_ptr->indices;
      if (*index_ptr < *min_occupied_D_ptr) {
        *min_occupied_D_ptr = *index_ptr;
      }
      if (*index_ptr > *max_occupied_D_ptr) {
        *max_occupied_D_ptr = *index_ptr;
      }
      ++index_ptr;
      if (*index_ptr < *min_occupied_H_ptr) {
        *min_occupied_H_ptr = *index_ptr;
      }
      if (*index_ptr > *max_occupied_H_ptr) {
        *max_occupied_H_ptr = *index_ptr;
      }
      ++index_ptr;
      if (*index_ptr < *min_occupied_P_ptr) {
        *min_occupied_P_ptr = *index_ptr;
      }
      if (*index_ptr > *max_occupied_P_ptr) {
        *max_occupied_P_ptr = *index_ptr;
      }
      ++index_ptr;
      if (*index_ptr < *min_occupied_M_ptr) {
        *min_occupied_M_ptr = *index_ptr;
      }
      if (*index_ptr > *max_occupied_M_ptr) {
        *max_occupied_M_ptr = *index_ptr;
      }
    }
    // Zero out index_of_potentially_occupied_mutant_class so we can use it
    // next time
    marray_uint_set(index_of_potentially_occupied_mutant_class,
                    class_indices_ptr->indices, 0);
  }
}

void compute_marginalized_driver_counts(marray_uint *counts, int Md, 
                                    int min_occupied_D, int min_occupied_H,
                                    int min_occupied_P, int min_occupied_M,
                                    int max_occupied_D, int max_occupied_H,
                                    int max_occupied_P, int max_occupied_M, 
                                    size_t *indices,
                                    marray_uint *marginalized_driver_counts, 
                                    size_t *marginalized_indices, FILE *fp) {
  unsigned int count_in_class, count;
  for (indices[0] = min_occupied_D; indices[0] <= max_occupied_D; indices[0]++) {
    marginalized_indices[1] = indices[0];
    count = 0;
    for (indices[1] = min_occupied_H; indices[1] <= max_occupied_H; indices[1]++) {
      for (indices[2] = min_occupied_P; indices[2] <= max_occupied_P; indices[2]++) {
        for (indices[3] = min_occupied_M; indices[3] <= max_occupied_M; indices[3]++) {
          count_in_class = marray_uint_get(counts, indices);
          count += count_in_class;
        }
      }
    }
    if (count > 0) {
      marray_uint_set(marginalized_driver_counts, marginalized_indices, count);
    }
  }
  fprintf(fp,"%d", (int)marginalized_indices[0]);
  for (marginalized_indices[1] = 0; marginalized_indices[1] <= Md;
      ++marginalized_indices[1]) {
    count = marray_uint_get(marginalized_driver_counts,
                            marginalized_indices);
      fprintf(fp, ",%u", count);
  }
  fprintf(fp,"\n");
  fflush(fp);
}

void compute_marginalized_housekeeper_counts(marray_uint *counts, int Mh, 
                                    int min_occupied_D, int min_occupied_H,
                                    int min_occupied_P, int min_occupied_M,
                                    int max_occupied_D, int max_occupied_H,
                                    int max_occupied_P, int max_occupied_M, 
                                    size_t *indices,
                                    marray_uint *marginalized_housekeeper_counts, 
                                    size_t *marginalized_indices, FILE *fp) {
  unsigned int count_in_class, count;
  for (indices[1] = min_occupied_H; indices[1] <= max_occupied_H; indices[1]++) {
    marginalized_indices[1] = indices[1];
    count = 0;
    for (indices[0] = min_occupied_D; indices[0] <= max_occupied_D; indices[0]++) {
      for (indices[2] = min_occupied_P; indices[2] <= max_occupied_P; indices[2]++) {
        for (indices[3] = min_occupied_M; indices[3] <= max_occupied_M; indices[3]++) {
          count_in_class = marray_uint_get(counts, indices);
          count += count_in_class;
        }
      }
    }
    if (count > 0) {
      marray_uint_set(marginalized_housekeeper_counts, marginalized_indices, count);
    }
  }
  fprintf(fp,"%d", (int)marginalized_indices[0]);
  for (marginalized_indices[1] = 0; marginalized_indices[1] <= Mh;
      ++marginalized_indices[1]) {
    count = marray_uint_get(marginalized_housekeeper_counts,
                            marginalized_indices);
      fprintf(fp, ",%u", count);
  }
  fprintf(fp,"\n");
  fflush(fp);
}

void compute_marginalized_passenger_counts(marray_uint *counts, int Mp, 
                                    int min_occupied_D, int min_occupied_H,
                                    int min_occupied_P, int min_occupied_M,
                                    int max_occupied_D, int max_occupied_H,
                                    int max_occupied_P, int max_occupied_M, 
                                    size_t *indices,
                                    marray_uint *marginalized_passenger_counts, 
                                    size_t *marginalized_indices, FILE *fp) {
  unsigned int count_in_class, count;
  for (indices[2] = min_occupied_P; indices[2] <= max_occupied_P; indices[2]++) {
    marginalized_indices[1] = indices[2];
    count = 0;
    for (indices[0] = min_occupied_D; indices[0] <= max_occupied_D; indices[0]++) {
      for (indices[1] = min_occupied_H; indices[1] <= max_occupied_H; indices[1]++) {
        for (indices[3] = min_occupied_M; indices[3] <= max_occupied_M; indices[3]++) {
          count_in_class = marray_uint_get(counts, indices);
          count += count_in_class;
        }
      }
    }
    if (count > 0) {
      marray_uint_set(marginalized_passenger_counts, marginalized_indices, count);
    }
  }
  fprintf(fp,"%d", (int)marginalized_indices[0]);
  for (marginalized_indices[1] = 0; marginalized_indices[1] <= Mp;
      ++marginalized_indices[1]) {
    count = marray_uint_get(marginalized_passenger_counts,
                            marginalized_indices);
      fprintf(fp, ",%u", count);
  }
  fprintf(fp,"\n");
  fflush(fp);
}

void compute_marginalized_mutator_counts(marray_uint *counts, int Mm, 
                                    int min_occupied_D, int min_occupied_H,
                                    int min_occupied_P, int min_occupied_M,
                                    int max_occupied_D, int max_occupied_H,
                                    int max_occupied_P, int max_occupied_M, 
                                    size_t *indices,
                                    marray_uint *marginalized_mutator_counts, 
                                    size_t *marginalized_indices, FILE *fp) {
  unsigned int count_in_class, count;
  for (indices[3] = min_occupied_M; indices[3] <= max_occupied_M; indices[3]++) {
    marginalized_indices[1] = indices[3];
    count = 0;
    for (indices[0] = min_occupied_D; indices[0] <= max_occupied_D; indices[0]++) {
      for (indices[1] = min_occupied_H; indices[1] <= max_occupied_H; indices[1]++) {
        for (indices[2] = min_occupied_P; indices[2] <= max_occupied_P; indices[2]++) {
          count_in_class = marray_uint_get(counts, indices);
          count += count_in_class;
        }
      }
    }
    if (count > 0) {
      marray_uint_set(marginalized_mutator_counts, marginalized_indices, count);
    }
  }
  fprintf(fp,"%d", (int)marginalized_indices[0]);
  for (marginalized_indices[1] = 0; marginalized_indices[1] <= Mm;
      ++marginalized_indices[1]) {
    count = marray_uint_get(marginalized_mutator_counts,
                            marginalized_indices);
      fprintf(fp, ",%u", count);
  }
  fprintf(fp,"\n");
  fflush(fp);
}

void compute_averages(marray_uint *counts, unsigned int N, 
                      int min_occupied_D, int min_occupied_H,
                      int min_occupied_P, int min_occupied_M,
                      int max_occupied_D, int max_occupied_H,
                      int max_occupied_P, int max_occupied_M, 
                      gsl_matrix *dataM, int index, size_t *indices) {
	
	// calc mean and variance of mutation levels in each class
	unsigned int meanM,meanD,meanP,meanH,varD,varH,varP,varM,count,tmp;
  double reciprocal_N, reciprocal_Nsquared;

	meanM = meanP = meanH = meanD = 0;
	varM = varP = varH = varD = 0;
  for (indices[0] = min_occupied_D; indices[0] <= max_occupied_D; indices[0]++) {
    for (indices[1] = min_occupied_H; indices[1] <= max_occupied_H; indices[1]++) {
      for (indices[2] = min_occupied_P; indices[2] <= max_occupied_P; indices[2]++) {
        for (indices[3] = min_occupied_M; indices[3] <= max_occupied_M; indices[3]++) {
          count = marray_uint_get(counts, indices);
          tmp = indices[0] * count;
          meanD += tmp;
          varD += indices[0] * tmp;
          tmp = indices[1] * count;
          meanH += tmp;
          varH += indices[1] * tmp;
          tmp = indices[2] * count;
          meanP += tmp;
          varP += indices[2] * tmp;
          tmp = indices[3] * count;
          meanM += tmp;
          varM += indices[3] * tmp;
        }
      }
    }
  }

  reciprocal_N = 1.0 / N;
  reciprocal_Nsquared = reciprocal_N * reciprocal_N;
	gsl_matrix_set(dataM, index, 0, gsl_matrix_get(dataM, index, 0) 
                                  + meanD * reciprocal_N);
	gsl_matrix_set(dataM, index, 1, gsl_matrix_get(dataM, index, 1) 
                                  + varD * reciprocal_N 
                                  - meanD*meanD * reciprocal_Nsquared);
	gsl_matrix_set(dataM, index, 2, gsl_matrix_get(dataM, index, 2) 
                                  + meanH * reciprocal_N);
	gsl_matrix_set(dataM, index, 3, gsl_matrix_get(dataM, index, 3) 
                                  + varH * reciprocal_N 
                                  - meanH*meanH * reciprocal_Nsquared);
	gsl_matrix_set(dataM, index, 4, gsl_matrix_get(dataM, index, 4) 
                                  + meanP * reciprocal_N);
	gsl_matrix_set(dataM, index, 5, gsl_matrix_get(dataM, index, 5) 
                                  + varP * reciprocal_N 
                                  - meanP*meanP * reciprocal_Nsquared);
	gsl_matrix_set(dataM, index, 6, gsl_matrix_get(dataM, index, 6) 
                                  + meanM * reciprocal_N);
	gsl_matrix_set(dataM, index, 7, gsl_matrix_get(dataM, index, 7) 
                                  + varM * reciprocal_N 
                                  - meanM*meanM * reciprocal_Nsquared);
	gsl_matrix_set(dataM, index, 8, gsl_matrix_get(dataM, index, 8) 
                                  + N);
	gsl_matrix_set(dataM, index, 9, gsl_matrix_get(dataM, index, 9) 
                                  + 1);
}

void report_data(int process_number, marray_uint *counts, unsigned int N, 
                  int rep, int gen, 
                  int min_occupied_D, int min_occupied_H, 
                  int min_occupied_P, int min_occupied_M,
                  int max_occupied_D, int max_occupied_H, 
                  int max_occupied_P, int max_occupied_M,
                  int time_first_sig_mut, int time_first_sig_ca, 
                  size_t *indices, DiversityIndices *diversity_indices_ptr,
                  FILE *fp) {
	
	int mostagg[4],mostmut[4],largest[4],least_domestic[4],most_hitchhiking[4],cancerP,mutP;
  int num_mostagg, num_mostmut, num_largest, num_least_domestic,
      num_most_hitchhiking;
	double prop_final_mut,prop_final_can;
	
	mostagg[0] = mostagg[1] = mostagg[2] = mostagg[3] = mostmut[0] = mostmut[1] = mostmut[2] = mostmut[3] = 0;
  largest[0] = largest[1] = largest[2] = largest[3] = 0;
  least_domestic[0] = least_domestic[1] = least_domestic[2] = least_domestic[3] 
    = 0;
  most_hitchhiking[0] = most_hitchhiking[1] = most_hitchhiking[2] = most_hitchhiking[3] = 0;
	
	// calc mean mutation levels
	double meanD,meanH,meanP,meanM,count,reciprocal_n;
  reciprocal_n = 1.0 / N;
	meanM = meanP = meanH = meanD = 0.0;
  for (indices[0] = min_occupied_D; indices[0] <= max_occupied_D; indices[0]++) {
    for (indices[1] = min_occupied_H; indices[1] <= max_occupied_H; indices[1]++) {
      for (indices[2] = min_occupied_P; indices[2] <= max_occupied_P; indices[2]++) {
        for (indices[3] = min_occupied_M; indices[3] <= max_occupied_M; indices[3]++) {
          count = (double)marray_uint_get(counts, indices) * reciprocal_n;
          meanD += indices[0] * count;
          meanH += indices[1] * count;
          meanP += indices[2] * count;
          meanM += indices[3] * count;
        }
      }
    }
  }
	
	num_mostagg = compute_most_aggressive_mutant_class(counts, 
                                      min_occupied_D, min_occupied_H,
                                      min_occupied_P, min_occupied_M, 
                                      max_occupied_D, max_occupied_H,
                                      max_occupied_P, max_occupied_M, 
                                      mostagg, indices);
	num_mostmut = compute_most_mutator_mutant_class(counts, 
                                    min_occupied_D, min_occupied_H,
                                    min_occupied_P, min_occupied_M, 
                                    max_occupied_D, max_occupied_H,
                                    max_occupied_P, max_occupied_M, 
                                    mostmut, indices);
	num_least_domestic = compute_least_domestic_mutant_class(counts, 
                                    min_occupied_D, min_occupied_H,
                                    min_occupied_P, min_occupied_M, 
                                    max_occupied_D, max_occupied_H,
                                    max_occupied_P, max_occupied_M, 
                                    least_domestic, indices);
	num_most_hitchhiking = compute_most_hitchhiking_mutant_class(counts, 
                                    min_occupied_D, min_occupied_H,
                                    min_occupied_P, min_occupied_M, 
                                    max_occupied_D, max_occupied_H,
                                    max_occupied_P, max_occupied_M, 
                                    most_hitchhiking, indices);
	num_largest = largest_mutant_class(counts, 
                      min_occupied_D, min_occupied_H, 
                      min_occupied_P, min_occupied_M, 
                      max_occupied_D, max_occupied_H, 
                      max_occupied_P, max_occupied_M, 
                      largest, indices);
	

	mutP = compute_number_of_mutator_cells(counts, 
                                          min_occupied_D, min_occupied_H, 
                                          min_occupied_P, min_occupied_M, 
                                          max_occupied_D, max_occupied_H, 
                                          max_occupied_P, max_occupied_M, 
                                          indices);
	cancerP = compute_number_of_cells_classed_as_cancerous(counts,
                min_occupied_D, min_occupied_H, min_occupied_P, min_occupied_M,
                max_occupied_D, max_occupied_H, max_occupied_P, max_occupied_M, 
                indices);
	
	prop_final_mut = (double) mutP * reciprocal_n;
	prop_final_can = (double) cancerP * reciprocal_n;

  compute_diversity_indices(counts,
                min_occupied_D, min_occupied_H, min_occupied_P, min_occupied_M,
                max_occupied_D, max_occupied_H, max_occupied_P, max_occupied_M, 
                indices, reciprocal_n, diversity_indices_ptr);
	
  fprintf(fp,"%d,%d,", process_number, ics.run_number);
	fprintf(fp,"%d,%d,%d,%d,", ics.gens,rep,gen,N);
  fprintf(fp,"%d,%d,%d,%d,", ics.Md,ics.Mh,ics.Mp,ics.Mm);
  fprintf(fp,"%d,%d,%d,%d,%d,%d,%d,%d,",
          min_occupied_D, max_occupied_D, 
          min_occupied_H, max_occupied_H,
          min_occupied_P, max_occupied_P,
          min_occupied_M, max_occupied_M);
  fprintf(fp,"%g,%g,%g,%g,%d,",
          ics.sD,ics.sH,ics.sM,ics.muB,ics.mostagg);
  fprintf(fp,"%g,%g,%g,%g,", meanD,meanH,meanP,meanM);
  fprintf(fp,"%d,%d,%d,%d,%d,",
          num_largest, largest[0],largest[1],largest[2],largest[3]);
  fprintf(fp,"%d,%d,%d,%d,%d,",
          num_mostagg, mostagg[0],mostagg[1],mostagg[2],mostagg[3]);
  fprintf(fp,"%d,%d,%d,%d,%d,",
          num_mostmut, mostmut[0],mostmut[1],mostmut[2],mostmut[3]);
  fprintf(fp,"%d,%d,%d,%d,%d,",
          num_least_domestic, least_domestic[0],least_domestic[1],least_domestic[2],least_domestic[3]);
  fprintf(fp,"%d,%d,%d,%d,%d,",
          num_most_hitchhiking, most_hitchhiking[0],most_hitchhiking[1],most_hitchhiking[2],most_hitchhiking[3]);
  fprintf(fp,"%d,%g,%d,%g,",
          time_first_sig_mut,prop_final_mut,time_first_sig_ca,prop_final_can);
  fprintf(fp,"%d,%g,%g\n",
          diversity_indices_ptr->richness,
          diversity_indices_ptr->shannon_index,
          diversity_indices_ptr->simpson_index);
	fflush(fp);
}

int compute_most_aggressive_mutant_class(marray_uint *counts, 
                                        int min_occupied_D, int min_occupied_H, 
                                        int min_occupied_P, int min_occupied_M,
                                        int max_occupied_D, int max_occupied_H, 
                                        int max_occupied_P, int max_occupied_M,
                                        int *class, size_t *indices) {
  unsigned int count;
  // We need signed loop indices so the loop condition i >= 0 will fail to be
  // satisfied at the end of the loop
  int i, j, k, l;
  for (i = max_occupied_D; i >= min_occupied_D; i--) {
    for (j = max_occupied_H; j >= min_occupied_H; j--) {
      for (k = max_occupied_P; k >= min_occupied_P; k--) {
        for (l = max_occupied_M; l >= min_occupied_M; l--) {
          indices[0] = i;
          indices[1] = j;
          indices[2] = k;
          indices[3] = l;
          count = marray_uint_get(counts, indices);
          if (count > 0) {
            // Since i is going down from max_occupied_D, this is the class with
            // the most driver mutations, i.e., the most aggressive class
            class[0] = indices[0];
            class[1] = indices[1];
            class[2] = indices[2];
            class[3] = indices[3];
            return (count);
          }
        }
      }
    }
  }
  return (0);
}

int compute_most_mutator_mutant_class(marray_uint *counts, 
                                      int min_occupied_D, int min_occupied_H, 
                                      int min_occupied_P, int min_occupied_M,
                                      int max_occupied_D, int max_occupied_H, 
                                      int max_occupied_P, int max_occupied_M,
                                      int *class, size_t *indices) {
  unsigned int count;
  // We need signed loop indices so the loop condition i >= 0 will fail to be
  // satisfied at the end of the loop
  int i, j, k, l;
  for (l = max_occupied_M; l >= min_occupied_M; l--) {
    for (i = max_occupied_D; i >= min_occupied_D; i--) {
      for (j = max_occupied_H; j >= min_occupied_H; j--) {
        for (k = max_occupied_P; k >= min_occupied_P; k--) {
          indices[0] = i;
          indices[1] = j;
          indices[2] = k;
          indices[3] = l;
          count = marray_uint_get(counts, indices);
          if (count > 0) {
            // Since indices[3] is going down from Mm, this is the class with
            // the most mutator mutations, i.e., the most mutator class
            class[0] = indices[0];
            class[1] = indices[1];
            class[2] = indices[2];
            class[3] = indices[3];
            return (count);
          }
        }
      }
    }
  }
  return (0);
}

int compute_least_domestic_mutant_class(marray_uint *counts, 
                                      int min_occupied_D, int min_occupied_H, 
                                      int min_occupied_P, int min_occupied_M,
                                      int max_occupied_D, int max_occupied_H, 
                                      int max_occupied_P, int max_occupied_M,
                                      int *class, size_t *indices) {
  unsigned int count;
  // We need signed loop indices so the loop condition i >= 0 will fail to be
  // satisfied at the end of the loop
  int i, j, k, l;
  for (j = max_occupied_H; j >= min_occupied_H; j--) {
    for (i = max_occupied_D; i >= min_occupied_D; i--) {
      for (k = max_occupied_P; k >= min_occupied_P; k--) {
        for (l = max_occupied_M; l >= min_occupied_M; l--) {
          indices[0] = i;
          indices[1] = j;
          indices[2] = k;
          indices[3] = l;
          count = marray_uint_get(counts, indices);
          if (count > 0) {
            // Since indices[1] is going down from Mh, this is the class with
            // the most housekeeper mutations, i.e., the least domestic class
            class[0] = indices[0];
            class[1] = indices[1];
            class[2] = indices[2];
            class[3] = indices[3];
            return (count);
          }
        }
      }
    }
  }
  return (0);
}

int compute_most_hitchhiking_mutant_class(marray_uint *counts, 
                                      int min_occupied_D, int min_occupied_H, 
                                      int min_occupied_P, int min_occupied_M,
                                      int max_occupied_D, int max_occupied_H, 
                                      int max_occupied_P, int max_occupied_M,
                                      int *class, size_t *indices) {
  unsigned int count;
  // We need signed loop indices so the loop condition i >= 0 will fail to be
  // satisfied at the end of the loop
  int i, j, k, l;
  for (k = max_occupied_P; k >= min_occupied_P; k--) {
    for (i = max_occupied_D; i >= min_occupied_D; i--) {
      for (j = max_occupied_H; j >= min_occupied_H; j--) {
        for (l = max_occupied_M; l >= min_occupied_M; l--) {
          indices[0] = i;
          indices[1] = j;
          indices[2] = k;
          indices[3] = l;
          count = marray_uint_get(counts, indices);
          if (count > 0) {
            // Since indices[2] is going down from Mp, this is the class with
            // the most passenger mutations, i.e., the most hitchhiking class
            class[0] = indices[0];
            class[1] = indices[1];
            class[2] = indices[2];
            class[3] = indices[3];
            return (count);
          }
        }
      }
    }
  }
  return (0);
}

int largest_mutant_class(marray_uint *counts, 
                        int min_occupied_D, int min_occupied_H, 
                        int min_occupied_P, int min_occupied_M,
                        int max_occupied_D, int max_occupied_H, 
                        int max_occupied_P, int max_occupied_M,
                        int *class, size_t *indices) {
  unsigned int count, max;
  max = 0;
  class[0] = class[1] = class[2] = class[3] = -1;
  for (indices[0] = min_occupied_D; indices[0] <= max_occupied_D; indices[0]++) {
    for (indices[1] = min_occupied_H; indices[1] <= max_occupied_H; indices[1]++) {
      for (indices[2] = min_occupied_P; indices[2] <= max_occupied_P; indices[2]++) {
        for (indices[3] = min_occupied_M; indices[3] <= max_occupied_M; indices[3]++) {
          count = marray_uint_get(counts, indices);
          if (count > max) {
            max = count;
            class[0] = indices[0];
            class[1] = indices[1];
            class[2] = indices[2];
            class[3] = indices[3];
          }
        }
      }
    }
  }
	return max;
}

int compute_number_of_mutator_cells(marray_uint *counts, 
                                    int min_occupied_D, int min_occupied_H, 
                                    int min_occupied_P, int min_occupied_M,
                                    int max_occupied_D, int max_occupied_H, 
                                    int max_occupied_P, int max_occupied_M,
                                    size_t *indices) {
	int count;
	count = 0;
	
	for(indices[0]=min_occupied_D; indices[0]<=max_occupied_D; indices[0]++) {
		for(indices[1]=min_occupied_H; indices[1]<=max_occupied_H; indices[1]++) {
      for(indices[2]=min_occupied_P; indices[2]<=max_occupied_P; indices[2]++) {
        for (indices[3]=GSL_MAX_INT(min_occupied_M,ics.mut_max);
              indices[3]<=max_occupied_M; indices[3]++) {
          count += marray_uint_get(counts, indices);
        }
      }
		}
	}
	return(count);
}

int compute_number_of_cells_classed_as_cancerous(marray_uint *counts,
                                                  int min_occupied_D, int min_occupied_H, 
                                                  int min_occupied_P, int min_occupied_M,
                                                  int max_occupied_D, int max_occupied_H, 
                                                  int max_occupied_P, int max_occupied_M,
                                                  size_t *indices) {
	int cells;
	cells = 0;
	for(indices[0]=GSL_MAX_INT(min_occupied_D, ics.mostagg); 
      indices[0]<=max_occupied_D; indices[0]++) {
		for(indices[1]=min_occupied_H; indices[1]<=max_occupied_H; indices[1]++) {
      for(indices[2]=GSL_MAX_INT(min_occupied_P,
                                ics.num_passengers_needed_for_cancer); 
          indices[2]<=max_occupied_P; indices[2]++) {
        for(indices[3]=min_occupied_M; indices[3]<=max_occupied_M; indices[3]++) {
          cells += marray_uint_get(counts, indices);
        }
      }
		}
	}
	return(cells);
}

void compute_diversity_indices(marray_uint *counts,
                              int min_occupied_D, int min_occupied_H, 
                              int min_occupied_P, int min_occupied_M,
                              int max_occupied_D, int max_occupied_H, 
                              int max_occupied_P, int max_occupied_M,
                              size_t *indices,
                              double reciprocal_population_size,
                              DiversityIndices *diversity_indices_ptr) {
  diversity_indices_ptr->richness = 0;
  diversity_indices_ptr->shannon_index = 0;
  diversity_indices_ptr->simpson_index = 0;
  unsigned int count;
  double proportion;
  for (indices[0] = min_occupied_D; indices[0] <= max_occupied_D; indices[0]++) {
    for (indices[1] = min_occupied_H; indices[1] <= max_occupied_H; indices[1]++) {
      for (indices[2] = min_occupied_P; indices[2] <= max_occupied_P; indices[2]++) {
        for (indices[3] = min_occupied_M; indices[3] <= max_occupied_M; indices[3]++) {
          count = marray_uint_get(counts, indices);
          if (count > 0) {
            ++(diversity_indices_ptr->richness);
            proportion = count * reciprocal_population_size;
            diversity_indices_ptr->shannon_index 
              -= proportion * log(proportion);
            diversity_indices_ptr->simpson_index += proportion * proportion;
          }
        }
      }
    }
  }
}
