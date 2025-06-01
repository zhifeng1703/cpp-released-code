#pragma once

#include <type_traits>
#include "colMat.hpp"

template <class ELEM_TYPE>
constexpr void skewl2m(ELEM_TYPE *MatA, INTE_TYPE lda, INTE_TYPE dim)
{
	for (auto ind = 0; ind < dim; MatA[ind * (lda + 1)] = 0, ind++)
		cblas_daxpby(dim - 1 - ind, -1.0, MatA + (ind * (lda + 1) + 1), 1, 0.0, MatA + (ind * (lda + 1) + lda), lda);
}

template <class ELEM_TYPE>
constexpr void skewl2m(ColMat<ELEM_TYPE> &MatA) { skewl2m(MatA.v, MatA.r, MatA.r); };

template <class ELEM_TYPE>
constexpr void skewl2m(View_ColMat<ELEM_TYPE> &ViewA) { skewl2m(ViewA.v, ViewA.ld, ViewA.r); };

template <class ELEM_TYPE>
constexpr ELEM_TYPE maxnrm(ELEM_TYPE *MatA, INTE_TYPE lda, INTE_TYPE row, INTE_TYPE col)
{
	ELEM_TYPE norm = 0;

	if (lda != row)
	{
		for (auto col_ind = 0; col_ind < col; col_ind++)
			for (auto row_ind = 0; row_ind < row; row_ind++)
				if (abs(MatA[row_ind + col_ind * lda]) > norm)
					norm = abs(MatA[row_ind + col_ind * lda]);
	}
	else
	{
		for (auto ind = 0; ind < col * row; ind++)
			if (abs(MatA[ind]) > norm)
				norm = abs(MatA[ind]);
	}
}

template <class ELEM_TYPE>
constexpr void maxnrm(ColMat<ELEM_TYPE> &MatA) { maxnrm(MatA.v, MatA.r, MatA.r, MatA.c); };

template <class ELEM_TYPE>
constexpr void maxnrm(View_ColMat<ELEM_TYPE> &ViewA) { maxnrm(ViewA.v, ViewA.ld, ViewA.r, ViewA.c); };