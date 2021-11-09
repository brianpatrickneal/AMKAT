/* Generate AMKAT observed statistics, with each using a separate application 
 * of AMKAT's filter method for feature selection
 
 AMKAT package for R
 Copyright (C) 2021, Brian Neal
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <RcppArmadillo.h>

#include "applyAmkatFilter.h"
#include "testSpearmanRho.h"
#include "getTailAreaSpearmanRho.h"
#include "computeSampleRanks.h"
#include "estimateSignalToNoise.h"
#include "generateKernelMatrix.h"

using namespace arma;

// 'x' and 'y' must have the same number of rows;
// lengths of 'y_variances' and of 'candidate_kernels' must both match 
// the column dimension of 'y';
// 'candidate_kernels' must contain values accepted by generateKernelMatrix;
// 'num_test_statistics' must be a strictly-positive integer
// see 'AMKAT/src/generateKernelMatrix.cpp';
// 'num_stats' must be a strictly-positive integer
// [[Rcpp::export]]
arma::vec generateTestStatMultiple(
    const arma::mat& y,
    const arma::vec& y_variances,
    const arma::mat& x,
    const Rcpp::CharacterVector& candidate_kernels,
    int num_test_statistics) {
  
  const int n = x.n_rows;
  const int num_kernels = candidate_kernels.size(); 
  const int num_y_variables = y.n_cols;
  arma::vec test_statistics(num_test_statistics, fill::zeros);
  arma::mat signal_to_noise(num_kernels, num_y_variables);
  arma::mat kernel_matrix(n, n);    
  uword index_of_max_snr;
  for (int k = 0; k < num_test_statistics; ++k) {
    arma::uvec selected_x_columns = applyAmkatFilter(y, x);
    for (int j = 0; j < num_kernels; ++j) {
      kernel_matrix = 
        generateKernelMatrix(x.cols(selected_x_columns), candidate_kernels[j]);
      for (int i = 0; i < num_y_variables; ++i) {
        signal_to_noise(j, i) = 
          estimateSignalToNoise(y.col(i), y_variances[i], kernel_matrix);
      }
    }
    for (int i = 0; i < num_y_variables; ++i) {
      index_of_max_snr = signal_to_noise.col(i).index_max();
      test_statistics[k] += signal_to_noise(index_of_max_snr, i);
    }
  }
  return test_statistics;
}