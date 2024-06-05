#pragma once

// This implementation computes the QR decomposition.

#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"

// LAPACKE_dorgqr(matrix_layout, m, m, p, a, lda, tau)

INTE_TYPE my_dorgqr(INTE_TYPE lapack_layout, INTE_TYPE m, INTE_TYPE n, INTE_TYPE k, REAL_TYPE *HHF, INTE_TYPE ldh, REAL_TYPE *tau);
INTE_TYPE my_dorgqr(INTE_TYPE m, INTE_TYPE n, INTE_TYPE k, REAL_TYPE *HHF, INTE_TYPE ldh, REAL_TYPE *tau);
INTE_TYPE my_dorgqr(INTE_TYPE m, INTE_TYPE n, INTE_TYPE k, REAL_TYPE *HHF, REAL_TYPE *tau);
INTE_TYPE my_dorgqr(INTE_TYPE k, ColMat<REAL_TYPE> &HHF, REAL_TYPE *tau);
INTE_TYPE my_dorgqr(INTE_TYPE k, View_ColMat<REAL_TYPE> &HHF, REAL_TYPE *tau);
INTE_TYPE my_dorgqr(ColMat<REAL_TYPE> &HHF, REAL_TYPE *tau);
INTE_TYPE my_dorgqr(View_ColMat<REAL_TYPE> &HHF, REAL_TYPE *tau);
