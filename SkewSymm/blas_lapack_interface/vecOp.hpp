#pragma once

#include <cassert>
#include <cstring>
#include "type_convention.hpp"

// #include "dnrm2.hpp"
#include "ddot.hpp"
#include "dscal.hpp"

//  dnrm2 is much slower than ddot as the former is designed for robust performance including pathological cases
//  see https://community.intel.com/t5/Intel-oneAPI-Math-Kernel-Library/cblas-dnrm2-much-slower-than-cblas-ddot/td-p/1022221

REAL_TYPE normsql2(INTE_TYPE n, const REAL_TYPE *v, INTE_TYPE incv);
REAL_TYPE normsql2(INTE_TYPE n, const REAL_TYPE *v);
REAL_TYPE norml2(INTE_TYPE n, const REAL_TYPE *v, INTE_TYPE incv);
REAL_TYPE norml2(INTE_TYPE n, const REAL_TYPE *v);

void scal(INTE_TYPE n, REAL_TYPE scale, REAL_TYPE *v);
