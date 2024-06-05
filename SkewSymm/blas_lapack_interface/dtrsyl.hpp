#pragma once

#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"
#include "matOp.hpp"

// This implemntation solves the Sylvester equation in form of AX+XB=C, where A, B are given upper quasi-triangular form from the real Schur decomposition.

INTE_TYPE my_dtrsyl(INTE_TYPE layout, CHAR_TYPE trana, CHAR_TYPE tranb, INTE_TYPE sgn, INTE_TYPE m, INTE_TYPE n,
                    REAL_TYPE *A, INTE_TYPE lda, REAL_TYPE *B, INTE_TYPE ldb, REAL_TYPE *C, INTE_TYPE ldc, REAL_TYPE *scale);

INTE_TYPE my_dtrsyl(CHAR_TYPE trana, CHAR_TYPE tranb, INTE_TYPE sgn, INTE_TYPE m, INTE_TYPE n,
                    REAL_TYPE *A, INTE_TYPE lda, REAL_TYPE *B, INTE_TYPE ldb, REAL_TYPE *C, INTE_TYPE ldc, REAL_TYPE *scale);

INTE_TYPE my_dtrsyl(CHAR_TYPE trana, CHAR_TYPE tranb, INTE_TYPE sgn, const ColMat<REAL_TYPE> &A, const ColMat<REAL_TYPE> &B, ColMat<REAL_TYPE> &C);
INTE_TYPE my_dtrsyl(CHAR_TYPE trana, CHAR_TYPE tranb, INTE_TYPE sgn, const View_ColMat<REAL_TYPE> &A, const View_ColMat<REAL_TYPE> &B, const View_ColMat<REAL_TYPE> &C);