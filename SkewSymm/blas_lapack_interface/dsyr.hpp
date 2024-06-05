#pragma once

#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"

// This implementation computes the symmetric rank-1 update A = alpha x * x' + A, where A is stored as ColMat.
// For the same computation but on the lower vector of A stored as LowTriMat with low_col_traversal, use dspr.

void my_dsyr(CBLAS_UPLO UPLO, INTE_TYPE n, REAL_TYPE alpha, REAL_TYPE *x, REAL_TYPE *a, INTE_TYPE lda);
void my_dsyr(INTE_TYPE n, REAL_TYPE alpha, REAL_TYPE *x, REAL_TYPE *a, INTE_TYPE lda);
void my_dsyr(INTE_TYPE n, REAL_TYPE alpha, REAL_TYPE *x, REAL_TYPE *a);
void my_dsyr(REAL_TYPE alpha, REAL_TYPE *x, const View_ColMat<REAL_TYPE>& ViewA);