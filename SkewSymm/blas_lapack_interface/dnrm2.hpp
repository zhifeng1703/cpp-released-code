#pragma once

// This implementation computes the l2 norm of a vector.

#include "mkl.h"
#include "type_convention.hpp"

REAL_TYPE my_dnrm2(INTE_TYPE n, const REAL_TYPE *v, INTE_TYPE incv);
REAL_TYPE my_dnrm2(INTE_TYPE n, const REAL_TYPE *v);