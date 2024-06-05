#pragma once
#include <cmath>
#include <iostream>
#include "mkl_types.h"
// #include "mkl.fi"

#define BOOL_TYPE int
#define REAL_TYPE double
#define INTE_TYPE MKL_INT
#define CMPX_TYPE MKL_Complex16
#define CHAR_TYPE char

#define LAYOUT_FLAG CBLAS_LAYOUT
#define TRANSP_FLAG CBLAS_TRANSPOSE

// CBLAS_TRANSPOSE transop(CBLAS_TRANSPOSE tran);

REAL_TYPE normofCMPX(const CMPX_TYPE &c);
void assignCMPX(CMPX_TYPE &c, const REAL_TYPE r, const REAL_TYPE i);
void assignCMPX(CMPX_TYPE &lhs, const CMPX_TYPE &rhs);
void additiCMPX(CMPX_TYPE &c, const CMPX_TYPE &a, const CMPX_TYPE &b);
void substrCMPX(CMPX_TYPE &c, const CMPX_TYPE &a, const CMPX_TYPE &b);
void multifCMPX(CMPX_TYPE &c, const CMPX_TYPE &a, const CMPX_TYPE &b);
void multifCMPX(CMPX_TYPE &c, REAL_TYPE a, const CMPX_TYPE &b);
void divideCMPX(CMPX_TYPE &c, const CMPX_TYPE &a, const CMPX_TYPE &b);
void inversCMPX(CMPX_TYPE &y, const CMPX_TYPE &x);
void exponeCMPX(CMPX_TYPE &y, const CMPX_TYPE &x);

std::ostream &operator<<(std::ostream &out, CMPX_TYPE c);