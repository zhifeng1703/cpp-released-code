#pragma once

#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"

// This implementation computes real vector y = alpha * A * x + beta * y
// using Intel(R) MKL function dgemm, where A is a matrix and , x and y 
// are vectors and alpha and beta are double precision scalars. 
// The matrix and vectors may be sent
// (1)  in the general array representation as follows:
//      (arr, row, col, leading_order, LAYOUT, TRANS);
// (2)  in the customized ColMat or View_ColMat in columnMajorMatrix.hpp

// Note that Users are respobible to manage the matrix size on their own,
// where A is in m x n, x is in n x 1 and y is in m x 1.

void my_dgemv(LAYOUT_FLAG LAYOUT, TRANSP_FLAG TRANS, INTE_TYPE m, INTE_TYPE n,
              REAL_TYPE alpha, REAL_TYPE *a_arr, INTE_TYPE lda, REAL_TYPE *x_arr, INTE_TYPE incx,
              REAL_TYPE beta, REAL_TYPE *y_arr, INTE_TYPE incy);

void my_dgemv(INTE_TYPE m, INTE_TYPE n, REAL_TYPE alpha, REAL_TYPE *a_arr, REAL_TYPE *x_arr, REAL_TYPE beta, REAL_TYPE *y_arr);

void my_dgemv(TRANSP_FLAG TRANS, REAL_TYPE alpha, const View_ColMat<REAL_TYPE> &A, const View_ColMat<REAL_TYPE> &ViewX,
              REAL_TYPE beta, const View_ColMat<REAL_TYPE> &ViewY);

void my_dgemv(TRANSP_FLAG TRANS, REAL_TYPE alpha, 
              const View_ColMat<REAL_TYPE> &A, const View_ColMat<REAL_TYPE> &ViewX, INTE_TYPE incx,
              REAL_TYPE beta, const View_ColMat<REAL_TYPE> &ViewY, INTE_TYPE incy);

void my_dgemv(REAL_TYPE alpha, const View_ColMat<REAL_TYPE> &A, const View_ColMat<REAL_TYPE> &ViewX,
              REAL_TYPE beta, const View_ColMat<REAL_TYPE> &ViewY);

void my_dgemv(REAL_TYPE alpha, const View_ColMat<REAL_TYPE> &A, const View_ColMat<REAL_TYPE> &ViewX, INTE_TYPE incx,
              REAL_TYPE beta, const View_ColMat<REAL_TYPE> &ViewY, INTE_TYPE incy);