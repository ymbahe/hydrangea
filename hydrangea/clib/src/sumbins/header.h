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

/**
 * @brief Structure to hold (pointers to) [NumPart] particle properties.
 */
struct st_part {
  long num;
  float* quant;
  int32_t* bin0;
  int32_t* bin1;
  int32_t* bin2;
  int32_t* bin3;
  int32_t* bin4;
  int32_t* bin5;
  int32_t* bin6;
  int32_t* bin7;
};

/**
  * @brief Structure to hold number of bins in each index dimension.
  */
struct st_bins {
  int32_t n0;
  int32_t n1;
  int32_t n2;
  int32_t n3;
  int32_t n4;
  int32_t n5;
  int32_t n6;
  int32_t n7;
};


/**
 * @brief Check that input is sensible and parse arguments if so.
 * 
 * If the input is inconsistent, a usage message is printed instead.
 * 
 * @param nArg Number of detected arguments.
 * @param argv[] The arguments passed to the main function.
 */
void get_input(int nArg, void* argv[], double** result);


static inline int64_t ind5d_bins(int64_t i, int64_t j, int64_t k, int64_t l,
				 int64_t m);

static inline int64_t ind6d_bins(int64_t i, int64_t j, int64_t k, int64_t l,
				 int64_t m, int64_t n);

static inline int64_t ind8d_bins(int64_t i, int64_t j, int64_t k, int64_t l,
				 int64_t m, int64_t n, int64_t o, int64_t p);


/* Global variables -- accessed everywhere, so easier to make global */

static struct st_part Part_struct;
static struct st_bins Bin_struct;

static struct st_part* part = &Part_struct;
static struct st_bins* bins = &Bin_struct;
