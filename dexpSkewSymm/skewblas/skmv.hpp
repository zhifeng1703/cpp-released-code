#pragma once
#include "skewMat.hpp"
#include "arrVec.hpp"

void skewblas_skmv(INTE_TYPE n, REAL_TYPE alpha, REAL_TYPE *A, INTE_TYPE lda, REAL_TYPE *x, INTE_TYPE incx, REAL_TYPE beta, REAL_TYPE *y, INTE_TYPE incy);
void skewblas_skmv(REAL_TYPE alpha, SkewSymmMat &A, REAL_TYPE *x, INTE_TYPE incx, REAL_TYPE beta, REAL_TYPE *y, INTE_TYPE incy);
void skewblas_skmv(REAL_TYPE alpha, SkewSymmMat &A, ArrVec<REAL_TYPE> &x, REAL_TYPE beta, ArrVec<REAL_TYPE> &y);
void skewblas_skmv(REAL_TYPE alpha, SkewSymmMat &A, View_ArrVec<REAL_TYPE> &x, REAL_TYPE beta, View_ArrVec<REAL_TYPE> &y);
void skewblas_skmv(REAL_TYPE alpha, ColMat<REAL_TYPE> &A, ArrVec<REAL_TYPE> &x, REAL_TYPE beta, ArrVec<REAL_TYPE> &y);
void skewblas_skmv(REAL_TYPE alpha, ColMat<REAL_TYPE> &A, View_ArrVec<REAL_TYPE> &x, REAL_TYPE beta, View_ArrVec<REAL_TYPE> &y);
