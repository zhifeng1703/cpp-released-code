#pragma once

// This implementation computes the pivoted-LU decomposition.


#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"

INTE_TYPE my_dgetrf(INTE_TYPE m, INTE_TYPE n, REAL_TYPE *A, INTE_TYPE lda, INTE_TYPE* P);
INTE_TYPE my_dgetrf(INTE_TYPE m, INTE_TYPE n, const View_ColMat<REAL_TYPE> &ViewA, INTE_TYPE* P);
