/* Estimates the standardized signal-to-noise ratio for variable y given its
 variance and a kernel similarity matrix for the independent variables

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

// Using mpf_float_100 type for increased precision (from boost header package)
// mpf_float types use mpfr and gmp libraries from BH (boost header) package;
// need flags for linking -lmpfr -lgmp to PKG_LIBS in makevars and makevars.win
#include <boost/multiprecision/mpfr.hpp>
namespace mp = boost::multiprecision;

using namespace arma;

// 'kernel_matrix' is a square matrix with the same row
// dimension as 'y', and that the length of 'y_variance' matches the column
// dimension of 'y'. Also assumes distribution of 'y' is centered around 0 or
// that the data has been centralized (e.g., by subtracting the sample mean from
// each value)
// [[Rcpp::export]]
double estimateSignalToNoise(const arma::vec& y,
                             double y_variance,
                             const arma::mat& kernel_matrix) {
  const int n = y.size();
  arma::mat kernel_matrix_diag0 = kernel_matrix;
  kernel_matrix_diag0.diag().zeros();
  const arma::mat centering_matrix = eye(n, n) - mat(n, n, fill::ones) / n;
  const arma::mat hk0 = centering_matrix * kernel_matrix_diag0;
  const arma::mat hk0hk0 = hk0 * hk0;
  const arma::mat hk0h = hk0 * centering_matrix;
  const arma::mat hk0h_hadamard = hk0h % hk0h;
  const mp::mpf_float_100 trace_hk0h_hadamard = trace(hk0h_hadamard);
  const mp::mpf_float_100 squared_trace_hk0 = trace(hk0) * trace(hk0);
  const mp::mpf_float_100 trace_hk0hk0 = trace(hk0hk0);
  const arma::vec y_standardized = y/sqrt(y_variance);
  const arma::vec fourth_power = pow(y_standardized, 4);
  const mp::mpf_float_100 n_float(n);
  const mp::mpf_float_100 fourth_moment = mean(fourth_power) - 3;
  const mp::mpf_float_100 snr_variance_float =
    (2 - 12 / (n_float - 1)) * trace_hk0hk0 -
    (2 / n_float) * squared_trace_hk0 +
    fourth_moment * ((6 / n_float) * trace_hk0hk0 +
    (1 / n_float) * squared_trace_hk0 + trace_hk0h_hadamard);
  const double snr_variance = snr_variance_float.convert_to<double>();
  const double signal_to_noise =
    (as_scalar(y.t() * kernel_matrix_diag0 * y)) / y_variance;
    return signal_to_noise / sqrt(snr_variance);
}