#pragma once

#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"

// This implementation computes the vector outer product A = alpha x * y' + A.

void my_dger(INTE_TYPE m, INTE_TYPE n, REAL_TYPE alpha, REAL_TYPE *x, REAL_TYPE *y, REAL_TYPE *a, INTE_TYPE lda);
void my_dger(INTE_TYPE m, INTE_TYPE n, REAL_TYPE alpha, REAL_TYPE *x, REAL_TYPE *y, const View_ColMat<REAL_TYPE> &ViewA);
