/* Generates the empirical centralized kernel similarity matrix for a set of 
 observations.  
 
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

using namespace Rcpp;

/* 'x' contains observations indexed by row.
 * 'kernel_function' accepts any one of the following strings for its value: 
 *   "gau", "lin", "quad", "exp", "IBS"
 * If adding to or modifying this list by changing the function definition, it is
 * also necessary to modify the variable 'programmed_kernels' defined in the 
 * main R function 'amkat' located in 'AMKAT/R/amkat_functions.R' */  
// [[Rcpp::export]]
arma::mat generateKernelMatrix(const arma::mat& x,
                               const Rcpp::String& kernel_function) {
  const arma::uword sample_size = x.n_rows;
  const int n(sample_size);
  const int p = x.n_cols;
  arma::mat kernel_matrix(n, n, arma::fill::zeros);
  if (kernel_function == "gau") {
    // Picking up gausskernel function from KRLS package
    Rcpp::Environment pkg = Environment::namespace_env("KRLS");
    Rcpp::Function f = pkg["gausskernel"];
    kernel_matrix = as<arma::mat>(f(x, p));
  } else if (kernel_function == "lin") {
    kernel_matrix = (x * x.t()) / p;
  } else if (kernel_function == "quad") {
    kernel_matrix = (x * x.t()) / p;
    kernel_matrix = pow(kernel_matrix + 1, 2.0);
  } else if (kernel_function == "exp") {
    for (arma::uword i = 0; i < sample_size; ++i) {
      const Rcpp::NumericVector x1 = 
        Rcpp::as<NumericVector>(Rcpp::wrap(x.row(i)));
      for (arma::uword j = 0; j < (i + 1); ++j) { // lower triangle
        const NumericVector x2 = Rcpp::as<NumericVector>(Rcpp::wrap(x.row(j)));
        kernel_matrix(i, j) = 
          exp((-sum(pow(x1, 2.0)) - 3 * sum(pow(x1 - x2, 2.0)) - 
          sum(pow(x2, 2.0))) / p);
      }
    }
    kernel_matrix = arma::symmatl(kernel_matrix); //reflect lower to upper
  } else if (kernel_function == "IBS") {
    for (arma::uword i = 0; i < sample_size; ++i) {
      const arma::rowvec x1 = x.row(i);
      for (arma::uword j = 0; j < (i + 1); ++j) { // lower triangle
        const arma::rowvec manhattan_dist_components(abs(x1 - x.row(j)));
        kernel_matrix(i, j) = 1 - (sum(manhattan_dist_components) / (2 * p));
      }
    }
    kernel_matrix = arma::symmatl(kernel_matrix); //reflect lower to upper
  }
  // empirical centralized kernel matrix
  arma::mat ker0 = kernel_matrix;
  ker0.diag().zeros();
  const arma::mat J(n, n, arma::fill::ones);
  return kernel_matrix - 
    (J * ker0 + ker0 * J - (J * ker0 * J) / n ) / (n - 1);
}