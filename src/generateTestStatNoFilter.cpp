/* Generate AMKAT observed test statistic with kernel selection results,
 * without using AMKAT's filter method for feature selection
 
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

#include "generateKernelMatrix.h"
#include "estimateSignalToNoise.h"

using namespace arma;

// NOTE: 'x' and 'y' must have the same number of rows;
// lengths of 'y_variances' and of 'candidate_kernels' must both match 
// the column dimension of 'y';
// 'candidate_kernels' must contain values accepted by generateKernelMatrix;
// see 'AMKAT/src/generateKernelMatrix.cpp'
// [[Rcpp::export]]
Rcpp::List generateTestStatNoFilter(
    const arma::mat& y,
    const arma::vec& y_variances,
    const arma::mat& x,
    const Rcpp::CharacterVector& candidate_kernels) {
  
  const int n = x.n_rows;
  const int num_kernels = candidate_kernels.size(); 
  const int num_y_variables = y.n_cols;        
  arma::mat kernel_matrix(n, n);
  arma::mat signal_to_noise(num_kernels, num_y_variables);
  uword index_of_max_snr;
  Rcpp::CharacterVector selected_kernels(num_y_variables);
  double test_statistic = 0;
  for (int j = 0; j < num_kernels; ++j) {
    kernel_matrix = generateKernelMatrix(x, candidate_kernels[j]);
    for (int i = 0; i < num_y_variables; ++i) {
      signal_to_noise(j, i) = 
        estimateSignalToNoise(y.col(i), y_variances[i], kernel_matrix);
    }
  }
  for (int i = 0; i < num_y_variables; ++i) {
    index_of_max_snr = signal_to_noise.col(i).index_max();
    selected_kernels[i] = candidate_kernels[index_of_max_snr];
    test_statistic += signal_to_noise(index_of_max_snr, i);
  }
  Rcpp::List output = 
    Rcpp::List::create(Rcpp::Named("test_statistic") = test_statistic,
                       Rcpp::Named("selected_kernels") = selected_kernels);
  return output;
}