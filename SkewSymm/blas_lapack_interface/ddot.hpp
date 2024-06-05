#pragma once

// This implementation computes the dot product x' * x

#include "mkl.h"
#include "type_convention.hpp"

REAL_TYPE my_ddot(INTE_TYPE n, const REAL_TYPE *x, INTE_TYPE incx, const REAL_TYPE *y, INTE_TYPE incy);
REAL_TYPE my_ddot(INTE_TYPE n, const REAL_TYPE *x, const REAL_TYPE *y);