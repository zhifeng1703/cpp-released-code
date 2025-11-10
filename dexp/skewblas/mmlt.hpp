#pragma once
#include "blasType.hpp"

INTE_TYPE _mmlt_blk_size(INTE_TYPE d);

void skewblas_mmlt(CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb, INTE_TYPE m, INTE_TYPE n, INTE_TYPE k,
                   REAL_TYPE alpha,
                   REAL_TYPE *A, INTE_TYPE lda,
                   REAL_TYPE *B, INTE_TYPE ldb,
                   REAL_TYPE beta,
                   REAL_TYPE *C, INTE_TYPE ldc);

void skewblas_mmlt(CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb, INTE_TYPE m, INTE_TYPE n, INTE_TYPE k,
                   REAL_TYPE alpha,
                   REAL_TYPE *A, INTE_TYPE lda,
                   REAL_TYPE *B, INTE_TYPE ldb,
                   REAL_TYPE beta,
                   REAL_TYPE *C, INTE_TYPE ldc, INTE_TYPE blk);