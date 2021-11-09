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


#ifndef AMKAT_SRC_GENERATETESTSTATMULTIPLE_H_
#define AMKAT_SRC_GENERATETESTSTATMULTIPLE_H_

arma::vec generateTestStatMultiple(
    const arma::mat& y,
    const arma::vec& y_variances,
    const arma::mat& x,
    const Rcpp::CharacterVector& candidate_kernels,
    int num_test_statistics) ;

#endif /* AMKAT_SRC_GENERATETESTSTATMULTIPLE_H_ */