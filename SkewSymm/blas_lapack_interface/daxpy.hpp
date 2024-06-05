#pragma once

// This implementation computes y = ax + y


#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"

void my_daxpy(INTE_TYPE n, REAL_TYPE a, REAL_TYPE *x, INTE_TYPE incx, REAL_TYPE *y, INTE_TYPE incy);
void my_daxpy(INTE_TYPE n, REAL_TYPE a, REAL_TYPE *x, REAL_TYPE *y);