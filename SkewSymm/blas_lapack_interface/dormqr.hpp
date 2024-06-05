#pragma once

// This implementation computes the matrix product C = QC, C = Q'C, C = CQ , C = CQ'
// where the orthogonal Q is characterized by a set of Househoulder reflectors (typically from dgeqef).

#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"

INTE_TYPE my_dormqr(INTE_TYPE lapack_layout, CHAR_TYPE lapack_side, CHAR_TYPE lapack_trans, INTE_TYPE m, INTE_TYPE n, INTE_TYPE k,
					const REAL_TYPE *HHF, INTE_TYPE ldh, const REAL_TYPE *tau, REAL_TYPE *C, INTE_TYPE ldc);
INTE_TYPE my_dormqr(CHAR_TYPE lapack_side, CHAR_TYPE lapack_trans, INTE_TYPE m, INTE_TYPE n, INTE_TYPE k,
					const REAL_TYPE *HHF, INTE_TYPE ldh, const REAL_TYPE *tau, REAL_TYPE *C, INTE_TYPE ldc);
INTE_TYPE my_dormqr(CHAR_TYPE lapack_side, CHAR_TYPE lapack_trans, const ColMat<REAL_TYPE> &HHF, const REAL_TYPE *tau, ColMat<REAL_TYPE> &C);
INTE_TYPE my_dormqr(CHAR_TYPE lapack_side, CHAR_TYPE lapack_trans, const View_ColMat<REAL_TYPE> &HHF, const REAL_TYPE *tau, View_ColMat<REAL_TYPE> &C);