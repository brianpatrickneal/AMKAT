/* Selects columns of 'x' that are related to one or more columns of 'y'

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

#include "testSpearmanRho.h"
#include "getTailAreaSpearmanRho.h"
#include "computeSampleRanks.h"

using namespace arma;

// NOTE: assumes 'y' and 'x' have the same number of rows
// [[Rcpp::export]]
arma::uvec applyAmkatFilter(const arma::mat& y,
                            const arma::mat& x) {
// for each column of x: tests Spearman's Rho with each column of y; if the
// minimum p-value is less than the value obtained using a permuted copy of x,
// the column is kept. If no columns of x are kept, the column with the lowest
// minimum p-value is kept by default
  const arma::uword sample_size = y.n_rows;
  const arma::uword num_x_variables = x.n_cols;
  const int n(sample_size);
  const int p(num_x_variables);
  const arma::mat x_permuted_rows(x.rows(randperm(n)));
  arma::vec min_pvalue(p, fill::ones);
  arma::vec min_pvalue_perm(p, fill::ones);
  double p_value;
  for (arma::uword i = 0; i < num_x_variables; ++i) {
    for (arma::uword j = 0; j < y.n_cols; ++j) {
      p_value = testSpearmanRho(y.col(j), x.col(i));
      min_pvalue[i] = std::min(p_value, min_pvalue[i]);
      p_value = testSpearmanRho(y.col(j), x_permuted_rows.col(i));
      min_pvalue_perm[i] = std::min(p_value, min_pvalue_perm[i]);
    }
  }
  arma::uvec selected_x_columns = find(min_pvalue < min_pvalue_perm);
  if (selected_x_columns.size() == 0) {
    arma::uvec default_x_column(1);
    default_x_column[1] = min_pvalue.index_min();
    return default_x_column;
  } else {
    return selected_x_columns;
  }
}