/* Computes the sample ranks of a numeric vector
 
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

using namespace arma;

// [[Rcpp::export]]
arma::vec computeSampleRanks(const arma::vec& x) {
  const arma::uword sample_size = x.size();
  const int n(sample_size);
  const arma::vec x_sorted = arma::sort(x);
  const arma::uvec x_sorted_indices = arma::sort_index(x);
  arma::vec x_ranks(n);
  arma::uword j = 0;
  for (arma::uword i = 0; i < sample_size; i = j + 1) {
    
    // how many values are tied with x_sorted[i]?
    j = i;
    while ((j < sample_size - 1) & (x_sorted[j] == x_sorted[j + 1])) { 
      j++; 
    }
    
    for (arma::uword k = i; k <= j; k++) {
      x_ranks[x_sorted_indices[k]] = (i + j + 2) / 2.; // averaging ties
    }
  }
  return x_ranks;
}