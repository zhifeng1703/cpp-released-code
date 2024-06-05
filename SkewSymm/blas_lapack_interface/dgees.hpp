#pragma once

#include "cmath"
#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"

// This implementation compute the real Schur decomposition of a real, non-symmetric matrix.

INTE_TYPE _dgees_select_nonzero(const REAL_TYPE *pr, const REAL_TYPE *pi);
INTE_TYPE _dgees_select_nonreal(const REAL_TYPE *pr, const REAL_TYPE *pi);

INTE_TYPE my_dgees(CHAR_TYPE JOBVS, CHAR_TYPE SORT, INTE_TYPE (*select)(const REAL_TYPE *, const REAL_TYPE *),
                   INTE_TYPE n, REAL_TYPE *a_arr, INTE_TYPE lda, INTE_TYPE *s_dim, REAL_TYPE *wr, REAL_TYPE *wi, REAL_TYPE *v_arr, INTE_TYPE ldv);