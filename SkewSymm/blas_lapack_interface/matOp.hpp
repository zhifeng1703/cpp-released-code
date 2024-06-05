#pragma once

#include <cassert>
#include <cstring>
#include <algorithm>
#include "type_convention.hpp"

#include "dnrm2.hpp"
#include "dscal.hpp"
#include "dgemm.hpp"
#include "dgemv.hpp"

void scal(INTE_TYPE n, INTE_TYPE k, REAL_TYPE s, REAL_TYPE *A, INTE_TYPE lda);
void scal(REAL_TYPE s, ColMat<REAL_TYPE> &A);
void scal(REAL_TYPE s, const View_ColMat<REAL_TYPE> &A);

void lmul(INTE_TYPE n, INTE_TYPE k, REAL_TYPE *A, INTE_TYPE lda, REAL_TYPE *B, INTE_TYPE ldb, REAL_TYPE *work);
void lmul(const ColMat<REAL_TYPE> &A, ColMat<REAL_TYPE> &B, REAL_TYPE *work);
// void lmul(REAL_TYPE a, ColMat<REAL_TYPE> &B);
// void lmul(const ColMat<REAL_TYPE> &A, REAL_TYPE* B, INTE_TYPE ldb, INTE_TYPE col);
// void lmul(const ColMat<REAL_TYPE> &A, REAL_TYPE* B, INTE_TYPE col);
// void lmul(const ColMat<REAL_TYPE> &A, REAL_TYPE* B);

void rmul(INTE_TYPE n, INTE_TYPE k, REAL_TYPE *A, INTE_TYPE lda, REAL_TYPE *B, INTE_TYPE ldb, REAL_TYPE *work);
void rmul(ColMat<REAL_TYPE> &A, const ColMat<REAL_TYPE> &B, REAL_TYPE *work);
// void rmul(ColMat<REAL_TYPE> &A, REAL_TYPE* B, INTE_TYPE ldb);
// void rmul(ColMat<REAL_TYPE> &A, REAL_TYPE* B);

void ladd(INTE_TYPE n, INTE_TYPE k, REAL_TYPE *A, INTE_TYPE lda, REAL_TYPE *B, INTE_TYPE ldb);
void ladd(const ColMat<REAL_TYPE> &A, ColMat<REAL_TYPE> &B);
// It is possible to further optimize the addtion and subtraction with calles BLAS level-1 operation but it seems unnecessary at this point.
// This implementation use the raw nested loops implementation and let complier optimizer to do their work.

void radd(INTE_TYPE n, INTE_TYPE k, REAL_TYPE *A, INTE_TYPE lda, REAL_TYPE *B, INTE_TYPE ldb);
void radd(ColMat<REAL_TYPE> &A, const ColMat<REAL_TYPE> &B);

void lsub(INTE_TYPE n, INTE_TYPE k, REAL_TYPE *A, INTE_TYPE lda, REAL_TYPE *B, INTE_TYPE ldb);
void lsub(const ColMat<REAL_TYPE> &A, ColMat<REAL_TYPE> &B);

void rsub(INTE_TYPE n, INTE_TYPE k, REAL_TYPE *A, INTE_TYPE lda, REAL_TYPE *B, INTE_TYPE ldb);
void rsub(ColMat<REAL_TYPE> &A, const ColMat<REAL_TYPE> &B);

REAL_TYPE normsqFro(const ColMat<REAL_TYPE> &A);
REAL_TYPE normsqFro(const View_ColMat<REAL_TYPE> &A);
REAL_TYPE normFro(const ColMat<REAL_TYPE> &A);
REAL_TYPE normFro(const View_ColMat<REAL_TYPE> &A);

// void axpby(INTE_TYPE m, INTE_TYPE n, INTE_TYPE k, REAL_TYPE alpha, REAL_TYPE* A, INTE_TYPE lda, REAL_TYPE beta, REAL_TYPE B, INTE_TYPE ldb, REAL_TYPE C, INTE_TYPE ldc, REAL_TYPE *work);
// void axpby(REAL_TYPE alpha, const ColMat<REAL_TYPE> &A, REAL_TYPE beta, const ColMat<REAL_TYPE> &B, ColMat<REAL_TYPE> &C, REAL_TYPE *work);

// congruence compute N = op(C) op(M) op(C') a necessary workspace W for op(C) op(M)
void congruence(CBLAS_TRANSPOSE tranc, CBLAS_TRANSPOSE tranm, INTE_TYPE m, INTE_TYPE n, REAL_TYPE *M, INTE_TYPE ldm, REAL_TYPE *C, INTE_TYPE ldc, REAL_TYPE *N, INTE_TYPE ldn, REAL_TYPE *W, INTE_TYPE ldw);
void congruence(CBLAS_TRANSPOSE tranc, CBLAS_TRANSPOSE tranm, const ColMat<REAL_TYPE> &M, const ColMat<REAL_TYPE> &C, ColMat<REAL_TYPE> &N, ColMat<REAL_TYPE> &W);
void congruence(CBLAS_TRANSPOSE tranc, CBLAS_TRANSPOSE tranm, const View_ColMat<REAL_TYPE> &M, const View_ColMat<REAL_TYPE> &C, const View_ColMat<REAL_TYPE> &N, const View_ColMat<REAL_TYPE> &W);
void congruence(CBLAS_TRANSPOSE tranc, INTE_TYPE m, INTE_TYPE n, REAL_TYPE *M, INTE_TYPE ldm, REAL_TYPE *C, INTE_TYPE ldc, REAL_TYPE *N, INTE_TYPE ldn, REAL_TYPE *W, INTE_TYPE ldw);
void congruence(CBLAS_TRANSPOSE tranc, const ColMat<REAL_TYPE> &M, const ColMat<REAL_TYPE> &C, ColMat<REAL_TYPE> &N, ColMat<REAL_TYPE> &W);
void congruence(CBLAS_TRANSPOSE tranc, const View_ColMat<REAL_TYPE> &M, const View_ColMat<REAL_TYPE> &C, const View_ColMat<REAL_TYPE> &N, const View_ColMat<REAL_TYPE> &W);
