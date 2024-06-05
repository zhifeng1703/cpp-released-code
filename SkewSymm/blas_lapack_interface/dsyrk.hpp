#pragma once

#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"

// This implementation computes the symmetric rank-k update A = alpha x * x' + beta A, where x is n x k and A is stored as ColMat
// For the same computation but on the lower vector of A stored as LowTriMat with low_col_traversal, use dsprk.

void my_dsyrk(CBLAS_UPLO UPLO, CBLAS_TRANSPOSE TRANS, INTE_TYPE n, INTE_TYPE k,
              REAL_TYPE alpha, REAL_TYPE *x, INTE_TYPE ldx, 
              REAL_TYPE beta, REAL_TYPE *a, INTE_TYPE lda);

void my_dsyrk(CBLAS_TRANSPOSE TRANS, INTE_TYPE n, INTE_TYPE k,
              REAL_TYPE alpha, REAL_TYPE *x, INTE_TYPE ldx, 
              REAL_TYPE beta, REAL_TYPE *a, INTE_TYPE lda);

void my_dsyrk(INTE_TYPE n, INTE_TYPE k,
              REAL_TYPE alpha, REAL_TYPE *x, INTE_TYPE ldx,
              REAL_TYPE beta, REAL_TYPE *a, INTE_TYPE lda);

void my_dsyrk(CBLAS_UPLO UPLO, INTE_TYPE n, INTE_TYPE k,
              REAL_TYPE alpha, REAL_TYPE *x, INTE_TYPE ldx,
              REAL_TYPE beta, REAL_TYPE *a, INTE_TYPE lda);

void my_dsyrk(INTE_TYPE n, INTE_TYPE k,
              REAL_TYPE alpha, const View_ColMat<REAL_TYPE> &ViewX,
              REAL_TYPE beta, const View_ColMat<REAL_TYPE> &ViewA);

void my_dsyrk(INTE_TYPE n, INTE_TYPE k,
              REAL_TYPE alpha, REAL_TYPE *x, INTE_TYPE ldx,
              REAL_TYPE beta, const View_ColMat<REAL_TYPE> &ViewA);

void my_dsyrk(REAL_TYPE alpha, REAL_TYPE *x, INTE_TYPE ldx, INTE_TYPE k,
              REAL_TYPE beta, const View_ColMat<REAL_TYPE> &ViewA);

void my_dsyrk(INTE_TYPE n, INTE_TYPE k,
              REAL_TYPE alpha, const View_ColMat<REAL_TYPE> &ViewX,
              REAL_TYPE beta, REAL_TYPE *a, INTE_TYPE lda);

void my_dsyrk(REAL_TYPE alpha, const View_ColMat<REAL_TYPE> &ViewX,
              REAL_TYPE beta, REAL_TYPE *a, INTE_TYPE lda);

void my_dsyrk(REAL_TYPE alpha, const View_ColMat<REAL_TYPE> &ViewX,
              REAL_TYPE beta, const View_ColMat<REAL_TYPE> &ViewA);
