/* Computes a left-tailed or right-tailed area for the sampling distribution 
 * of Spearman's rho
 
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

#ifdef __cplusplus
extern "C" {
#endif
#include "computeTailAreaSpearmanRho.h"
#ifdef __cplusplus
}
#endif
#include <Rcpp.h>

// [[Rcpp::export]]
double getTailAreaSpearmanRho(double q, int n, int ltail) {
  const double ifault = 0;
  double pval;
  computeTailAreaSpearmanRho(n, q, ifault, ltail, &pval);
  return pval;
}