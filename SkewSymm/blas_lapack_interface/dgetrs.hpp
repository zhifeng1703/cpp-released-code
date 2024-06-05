#pragma once

// This implementation computes the pivoted-LU decomposition.


#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"

// lapack_int LAPACKE_dgetrs (int matrix_layout , char trans , lapack_int n , lapack_int
// nrhs , const double * a , lapack_int lda , const lapack_int * ipiv , double * b ,
// lapack_int ldb )

INTE_TYPE my_dgetrs(char TRANSA, INTE_TYPE d, INTE_TYPE c, REAL_TYPE *LU, INTE_TYPE ldlu, INTE_TYPE *P, REAL_TYPE *X, INTE_TYPE ldx);
INTE_TYPE my_dgetrs(char TRANSA, INTE_TYPE d, REAL_TYPE *LU, INTE_TYPE ldlu, INTE_TYPE *P, const View_ColMat<REAL_TYPE> &ViewX);
INTE_TYPE my_dgetrs(char TRANSA, INTE_TYPE *P, const View_ColMat<REAL_TYPE> &ViewLU, const View_ColMat<REAL_TYPE> &ViewX);