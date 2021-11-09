/* Modified version of r-source/src/library/stats/src/prho.c
  
 * Computes a left-tailed or right-tailed area for the sampling distribution
 *  of Spearman's rho.
 
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2000-2016     The R Core Team
 *  Copyright (C) 2003		The R Foundation
 *  based on AS 89 (C) 1975 Royal Statistical Society
 
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
 
 * Modifications by Brian Neal [January, 2021]:
 *  - function name changed from 'prho' to 'computePvalueSpearmanRho'
 *  - output parameter 'pv' now placed after input parameters in argument order
 *  - parameter 'double *pv' reformatted as 'double* pv' in the function 
 *    declaration to emphasize the use of * for pointer declaration and avoid 
 *    confusion with the dereference operator *, since *pv is used within the 
 *    function definition to dereference the value pointed to by pv
 *  - removed the second function, 'pRho', used to call prho from R  
 */

#ifndef AMKAT_SRC_COMPUTETAILAREASPEARMANRHO_H_
#define AMKAT_SRC_COMPUTETAILAREASPEARMANRHO_H_

void computeTailAreaSpearmanRho(int n, double is, int ifault, int lower_tail, double* pv);

#endif /* AMKAT_SRC_COMPUTETAILAREASPEARMANRHO_H_ */