/*
    Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright © 2011-2017, Xiang Zhou
    Copyright © 2017, Peter Carbonetto
    Copyright © 2017, Pjotr Prins

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __FASTOPENBLAS_H__
#define __FASTOPENBLAS_H__

#include <assert.h>
#include <iostream>

// Include GSL headers first (they define CBLAS enums as regular enums)
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_cblas.h"

// Don't include system cblas.h as it conflicts with GSL's cblas definitions
// extern "C"
// {
//    #include <cblas.h>   // For OpenBlas / Atlas
// }

void fast_cblas_dgemm(const enum CBLAS_ORDER Order,
                      const enum CBLAS_TRANSPOSE TransA,
                      const enum CBLAS_TRANSPOSE TransB,
                      const size_t M,
                      const size_t N,
                      const size_t K,
                      const double alpha,
                      const double *A,
                      const size_t lda,
                      const double *B,
                      const size_t ldb,
                      const double beta,
                      double *C,
                      const size_t ldc);

#endif // __FASTOPENBLAS_H_
