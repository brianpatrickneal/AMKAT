/* Computes the p-value for a a two-tailed test of H_0: Spearman's Rho = 0

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

#include "getTailAreaSpearmanRho.h"
#include "computeSampleRanks.h"

using namespace Rcpp;

// The implementation is based on the one found in
// r-source/src/library/stats/R/cor.test.R
// NOTE: assumes 'x' and 'y' have the same length
// [[Rcpp::export]]
double testSpearmanRho(const arma::vec& x,
                       const arma::vec& y) {
  const int n = x.n_elem;
  const arma::colvec x_ranks = computeSampleRanks(x);
  const arma::colvec y_ranks = computeSampleRanks(y);
  const double rho_spearman = arma::as_scalar(arma::cor(x_ranks, y_ranks));
  const int left_tailed = (rho_spearman > 0);

  // Are there ties in either of x or y?
  Rcpp::NumericVector x_unique_values = Rcpp::wrap(x);
  Rcpp::NumericVector y_unique_values = Rcpp::wrap(y);
  x_unique_values = Rcpp::unique(x_unique_values);
  y_unique_values = Rcpp::unique(y_unique_values);
  const double x_num_unique_values = x_unique_values.size();
  const double y_num_unique_values = y_unique_values.size();
  const bool ties = (std::min(x_num_unique_values, y_num_unique_values) < n);
  double tail_area = 1;
  if ((n <= 1290) & !ties) {
    // using algorithm AS 89
    const double q = (n * (n * n - 1)) * (1 - rho_spearman) / 6;
    tail_area = getTailAreaSpearmanRho(round(q) + (2 * left_tailed),
                                       n, left_tailed);
  } else {
    // using Student's t
    // Rf_pt calls pt() in R, but is missing the ncp argument;
    // thus args are Rf_pt(quantile, df, lower.tail, log.p)
    tail_area =
      Rf_pt(rho_spearman / sqrt((1 - rho_spearman * rho_spearman) / (n - 2)),
            n - 2, !left_tailed, 0);
  }
  return (std::min(1.0, 2 * tail_area));
}