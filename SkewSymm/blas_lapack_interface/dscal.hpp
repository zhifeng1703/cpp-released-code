#pragma once

// This implementation scales the vector with a constant

#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"

void my_dscal(INTE_TYPE n, REAL_TYPE s, REAL_TYPE *x, INTE_TYPE incx);
void my_dscal(INTE_TYPE n, REAL_TYPE s, REAL_TYPE *x);