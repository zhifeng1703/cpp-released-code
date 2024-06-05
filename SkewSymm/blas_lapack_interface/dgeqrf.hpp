#pragma once

// This implementation computes the QR decomposition.

#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"

INTE_TYPE my_dgeqrf(INTE_TYPE lapack_layout, INTE_TYPE m, INTE_TYPE n, REAL_TYPE *A, INTE_TYPE lda, REAL_TYPE *tau);
INTE_TYPE my_dgeqrf(INTE_TYPE m, INTE_TYPE n, REAL_TYPE *A, INTE_TYPE lda, REAL_TYPE *tau);
INTE_TYPE my_dgeqrf(INTE_TYPE m, INTE_TYPE n, REAL_TYPE *A, REAL_TYPE *tau);
INTE_TYPE my_dgeqrf(ColMat<REAL_TYPE> &A, REAL_TYPE *tau);
INTE_TYPE my_dgeqrf(View_ColMat<REAL_TYPE> &A, REAL_TYPE *tau);
