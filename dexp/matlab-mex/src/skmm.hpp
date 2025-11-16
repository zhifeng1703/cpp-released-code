#pragma once
#include "skewMat.hpp"
#include "arrVec.hpp"

void skewblas_skmm(CHAR_TYPE side, INTE_TYPE m, INTE_TYPE n, REAL_TYPE alpha, REAL_TYPE *A, INTE_TYPE lda, REAL_TYPE *B, INTE_TYPE ldb, REAL_TYPE beta, REAL_TYPE *C, INTE_TYPE ldc);

void skewblas_skmm(CHAR_TYPE side, REAL_TYPE alpha, SkewSymmMat &A, ColMat<REAL_TYPE> &B, REAL_TYPE beta, ColMat<REAL_TYPE> &C);
void skewblas_skmm(CHAR_TYPE side, REAL_TYPE alpha, SkewSymmMat &A, View_ColMat<REAL_TYPE> &B, REAL_TYPE beta, ColMat<REAL_TYPE> &C);
void skewblas_skmm(CHAR_TYPE side, REAL_TYPE alpha, SkewSymmMat &A, ColMat<REAL_TYPE> &B, REAL_TYPE beta, View_ColMat<REAL_TYPE> &C);
void skewblas_skmm(CHAR_TYPE side, REAL_TYPE alpha, SkewSymmMat &A, View_ColMat<REAL_TYPE> &B, REAL_TYPE beta, View_ColMat<REAL_TYPE> &C);

void skewblas_skmm(CHAR_TYPE side, REAL_TYPE alpha, ColMat<REAL_TYPE> &A, ColMat<REAL_TYPE> &B, REAL_TYPE beta, ColMat<REAL_TYPE> &C);
void skewblas_skmm(CHAR_TYPE side, REAL_TYPE alpha, View_ColMat<REAL_TYPE> &A, ColMat<REAL_TYPE> &B, REAL_TYPE beta, ColMat<REAL_TYPE> &C);
void skewblas_skmm(CHAR_TYPE side, REAL_TYPE alpha, ColMat<REAL_TYPE> &A, View_ColMat<REAL_TYPE> &B, REAL_TYPE beta, ColMat<REAL_TYPE> &C);
void skewblas_skmm(CHAR_TYPE side, REAL_TYPE alpha, ColMat<REAL_TYPE> &A, ColMat<REAL_TYPE> &B, REAL_TYPE beta, View_ColMat<REAL_TYPE> &C);
void skewblas_skmm(CHAR_TYPE side, REAL_TYPE alpha, View_ColMat<REAL_TYPE> &A, View_ColMat<REAL_TYPE> &B, REAL_TYPE beta, ColMat<REAL_TYPE> &C);
void skewblas_skmm(CHAR_TYPE side, REAL_TYPE alpha, View_ColMat<REAL_TYPE> &A, ColMat<REAL_TYPE> &B, REAL_TYPE beta, View_ColMat<REAL_TYPE> &C);
void skewblas_skmm(CHAR_TYPE side, REAL_TYPE alpha, ColMat<REAL_TYPE> &A, View_ColMat<REAL_TYPE> &B, REAL_TYPE beta, View_ColMat<REAL_TYPE> &C);
void skewblas_skmm(CHAR_TYPE side, REAL_TYPE alpha, View_ColMat<REAL_TYPE> &A, View_ColMat<REAL_TYPE> &B, REAL_TYPE beta, View_ColMat<REAL_TYPE> &C);
