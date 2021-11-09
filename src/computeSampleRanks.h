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

#ifndef AMKAT_SRC_COMPUTESAMPLERANKS_H_
#define AMKAT_SRC_COMPUTESAMPLERANKS_H_

arma::vec computeSampleRanks(const arma::vec& x);

#endif /* AMKAT_SRC_COMPUTESAMPLERANKS_H_ */