#pragma once

#include "type_convention.hpp"
#include "dgeqrf.hpp"
#include "dorgqr.hpp"
#include "dormqr.hpp"

class HouseholderQrFactor
{
public:
	ColMat<REAL_TYPE> H;
	REAL_TYPE *t;
	INTE_TYPE r;
	INTE_TYPE c;

	HouseholderQrFactor() : H(ColMat<REAL_TYPE>()), t(nullptr), r(0), c(0){};
	HouseholderQrFactor(INTE_TYPE row, INTE_TYPE col) : H(ColMat<REAL_TYPE>(row, col)), t(new REAL_TYPE[col]), r(row), c(col){};
	void copy(const HouseholderQrFactor &src)
	{
		this->H.copy(src.H);
		memcpy(this->t, src.t, sizeof(REAL_TYPE) * c);
	};
	HouseholderQrFactor(const HouseholderQrFactor &src) : HouseholderQrFactor(src.r, src.c) { this->copy(src); };
	void swap(HouseholderQrFactor &rhs)
	{
		this->H.swap(rhs.H);
		using std::swap;
		swap(this->t, rhs.t);
		swap(this->r, rhs.r);
		swap(this->c, rhs.c);
	};
	HouseholderQrFactor &operator=(const HouseholderQrFactor &src)
	{
		HouseholderQrFactor temp(src);
		swap(temp);
		return (*this);
	}
	~HouseholderQrFactor()
	{
		if (t)
			delete[] t;
		t = nullptr;
	}

	void QR(const REAL_TYPE *A, INTE_TYPE lda)
	{
		H.assign(A, lda);
		my_dgeqrf(H, t);
	};
	void QR(const REAL_TYPE *A)
	{
		H.assign(A);
		my_dgeqrf(H, t);
	};
	void QR(const ColMat<REAL_TYPE> &A)
	{
		H.assign(A);
		my_dgeqrf(H, t);
	}
	void QR(const View_ColMat<REAL_TYPE> &A)
	{
		H.assign(A);
		my_dgeqrf(H, t);
	};

	void R(ColMat<REAL_TYPE> &MatR)
	{
		auto src_ptr = H.v;
		auto des_ptr = MatR.v;

		for (INTE_TYPE col_ind = 0; col_ind < c; col_ind++, src_ptr += H.r, des_ptr += MatR.r)
			memcpy(des_ptr, src_ptr, sizeof(REAL_TYPE) * (col_ind + 1));
	}
	void R(View_ColMat<REAL_TYPE> &MatR)
	{
		auto src_ptr = H.v;
		auto des_ptr = MatR.v;

		for (INTE_TYPE col_ind = 0; col_ind < c; col_ind++, src_ptr += H.r, des_ptr += MatR.ld)
			memcpy(des_ptr, src_ptr, sizeof(REAL_TYPE) * (col_ind + 1));
	};

	void Q(ColMat<REAL_TYPE> &MatQ, INTE_TYPE n, INTE_TYPE k)
	{
		// k is the number of Householder reflectors that forms Q_k
		// n is the number of the leading columns in Q_k being computed.
		// When n is omitted, the number of columns in A is used.
		// When k is omitted, the number of columns in H is used.
		memcpy(MatQ.v, H.v, sizeof(REAL_TYPE) * r * n);
		my_dorgqr(r, n, k, MatQ.v, t);
	};
	void Q(ColMat<REAL_TYPE> &MatQ, INTE_TYPE k)
	{
		memcpy(MatQ.v, H.v, sizeof(REAL_TYPE) * r * MatQ.c);
		my_dorgqr(k, MatQ, t);
	};
	void Q(ColMat<REAL_TYPE> &MatQ)
	{
		memcpy(MatQ.v, H.v, sizeof(REAL_TYPE) * r * c);
		my_dorgqr(c, MatQ, t);
	};
};
