/* -----------------------------------------------------------------------
# This file is part of the hydrangea tools package
# Copyright (C) 2019 Yannick Bahe (bahe@strw.leidenuniv.nl)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "assert.h"
#include "header.h"   /* Declares globals part and bin */

/**
 * @brief Main public function to bin a quantity by labels. 
 */
int sumbins(int argc, void *argv[]) {

  int64_t ii;
  double* result = NULL;
  double* kahan_c = NULL;

  /* Parse input parameters into internally-useful form */
  get_input(argc, argv, &result, &kahan_c);

  /* Go through particles and increment the right bin */
  for (ii = 0; ii < part->num; ii++) {
    const int64_t ind1d = ind8d_bins(part->bin0[ii], part->bin1[ii],
				     part->bin2[ii], part->bin3[ii],
				     part->bin4[ii], part->bin5[ii],
				     part->bin6[ii], part->bin7[ii]);
    const double kahan_y = (double) part->quant[ii] - kahan_c[ind1d];
    const double kahan_t = (double) result[ind1d] + kahan_y;
    kahan_c[ind1d] = (kahan_t - result[ind1d]) - kahan_y;
    result[ind1d] = kahan_t;
  }
  return 0;
}


/**
 * @brief Check that input is sensible and parse arguments if so.
 * 
 * If the input is inconsistent, a usage message is printed instead.
 * 
 * @param nArg Number of detected arguments.
 * @param argv[] The arguments passed to the main function.
 */
void get_input(int nArg, void* argv[], double** result, double** kahan_c) {
  
  /* Check that input is as expected and abort with instructions if not */
  if(nArg != 20) {
    printf("\n\nWrong number of arguments supplied to sum_index.\n\n");
    printf("Required arguments (19 total):\n"
	   "-------------------------------------------------------\n"
	   " 0     NumPart     (int64)  --> Number of particles    \n"
	   " 1-8   N_Bins[0-7] (int32)  --> Number of bins [0-7]   \n"
	   " 9     Quantity    (flt32)  --> Quantity to sum        \n"
	   " 10-17 Index[0-7]  (int32)  --> Bin index [0-7]        \n"
	   " 18    Result      (flt64)  --> Result array           \n"
	   " 19    Temporary   (flt64)  --> Temporary              \n"
	   "-------------------------------------------------------\n"
	   "\n\n");
    exit(0);
  }
  
  /* Number of particles in input data */
  part->num = *(long*) argv[0]; 

  bins->n0 = *(int32_t*) argv[1];
  bins->n1 = *(int32_t*) argv[2];
  bins->n2 = *(int32_t*) argv[3];
  bins->n3 = *(int32_t*) argv[4];
  bins->n4 = *(int32_t*) argv[5];
  bins->n5 = *(int32_t*) argv[6];
  bins->n6 = *(int32_t*) argv[7];
  bins->n7 = *(int32_t*) argv[8];

  part->quant = (float*) argv[9];

  part->bin0 = (int32_t*) argv[10];
  part->bin1 = (int32_t*) argv[11];
  part->bin2 = (int32_t*) argv[12];
  part->bin3 = (int32_t*) argv[13];
  part->bin4 = (int32_t*) argv[14];
  part->bin5 = (int32_t*) argv[15];
  part->bin6 = (int32_t*) argv[16];
  part->bin7 = (int32_t*) argv[17];

  *result = (double*) argv[18];
  *kahan_c = (double*) argv[19];
 
  return;
}


static inline int64_t ind5d_bins(int64_t i, int64_t j, int64_t k, int64_t l,
				 int64_t m) {

  const int64_t n1 = (int64_t) (bins->n0);
  const int64_t n2 = (int64_t) (bins->n1);
  const int64_t n3 = (int64_t) (bins->n2);
  const int64_t n4 = (int64_t) (bins->n3);
  const int64_t n5 = (int64_t) (bins->n4);

  assert(i < n1);
  assert(j < n2);
  assert(k < n3);
  assert(l < n4);
  assert(m < n5);

  return (n2*n3*n4*n5*i + n3*n4*n5*j + n4*n5*k + n5*l + m);
}

static inline int64_t ind6d_bins(int64_t i, int64_t j, int64_t k, 
				 int64_t l, int64_t m, int64_t n) {

  const int64_t n1 = (int64_t) (bins->n0);
  const int64_t n2 = (int64_t) (bins->n1);
  const int64_t n3 = (int64_t) (bins->n2);
  const int64_t n4 = (int64_t) (bins->n3);
  const int64_t n5 = (int64_t) (bins->n4);
  const int64_t n6 = (int64_t) (bins->n5);

  assert(i < n1);
  assert(j < n2);
  assert(k < n3);
  assert(l < n4);
  assert(m < n5);
  assert(n < n6);

  return (n2*n3*n4*n5*n6*i + n3*n4*n5*n6*j + n4*n5*n6*k + n5*n6*l + n6*m + n);
}

static inline int64_t ind8d_bins(int64_t i, int64_t j, int64_t k, 
				 int64_t l, int64_t m, int64_t n,
				 int64_t o, int64_t p) {

  const int64_t n1 = (int64_t) (bins->n0);
  const int64_t n2 = (int64_t) (bins->n1);
  const int64_t n3 = (int64_t) (bins->n2);
  const int64_t n4 = (int64_t) (bins->n3);
  const int64_t n5 = (int64_t) (bins->n4);
  const int64_t n6 = (int64_t) (bins->n5);
  const int64_t n7 = (int64_t) (bins->n6);
  const int64_t n8 = (int64_t) (bins->n7);

  assert(i < n1);
  assert(j < n2);
  assert(k < n3);
  assert(l < n4);
  assert(m < n5);
  assert(n < n6);
  assert(o < n7);
  assert(p < n8);

  return (n2*n3*n4*n5*n6*n7*n8*i + 
	  n3*n4*n5*n6*n7*n8*j + 
	  n4*n5*n6*n7*n8*k + 
	  n5*n6*n7*n8*l + n6*n7*n8*m + n7*n8*n + n8*o + p);
}


