#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

#include "mkl.h"
#include "blas_lapack_interface/type_convention.hpp"
#include "blas_lapack_interface/matOp.hpp"

#include "blas_lapack_interface/colMajMat.hpp"
#include "blas_lapack_interface/skewSymmMat.hpp"
#include "blas_lapack_interface/specOrthMat.hpp"

#include "blas_lapack_interface/expmpa.hpp"
#include "blas_lapack_interface/gmres.hpp"

#include "dexpSkewSymm.hpp"

#define _StGeoNewton_ExpZ_Scale_Order 4

struct StiefelGeoNewtonOperator
{
public:
	INTE_TYPE n;
	INTE_TYPE p;
	dexpSkewSymmPara *dexp; // dexp_S
	SkewSymmMat *Y;			// n x n Skew symmetric matrices Y = (dexp_S^{-1})[X]. It must carry blk_traversal
	SkewSymmMat *X;			// n x n Skew symmetric matrices in forms of S_{0, 0, Z}. It must carry a traversal.

	void Action(REAL_TYPE *y, REAL_TYPE *x) const
	{
		auto X_col = X->col(p);
		auto Y_col = Y->col(p);

		INTE_TYPE vec_ind = 0;
		for (INTE_TYPE col_ind = p; col_ind < n; col_ind++, X_col += n)
			for (INTE_TYPE row_ind = col_ind + 1; row_ind < n; row_ind++)
				X_col[row_ind] = x[vec_ind++];
		X->mat2vec();
		X->vec2mat();

		// X->printf("MatX:\n");
		dexpSkewSymm_inverse(*Y, *X, *dexp);
		// Y->printf("MatY:\n");

		vec_ind = 0;
		for (INTE_TYPE col_ind = p; col_ind < n; col_ind++, Y_col += n)
			for (INTE_TYPE row_ind = col_ind + 1; row_ind < n; row_ind++)
				y[vec_ind++] = Y_col[row_ind];

		// vec_ind = 0;
		// Y_col = Y->col(p);
		// for (INTE_TYPE col_ind = p; col_ind < n; col_ind++, Y_col += n)
		//	for (INTE_TYPE row_ind = col_ind + 1; row_ind < n; row_ind++)
		//		printf("%1.3f\t", Y_col[row_ind]);
		// printf("\n");
	};
	void Remainder(REAL_TYPE *r, REAL_TYPE *x) const
	{
		auto X_col = X->col(p);
		auto Y_col = Y->col(p);

		INTE_TYPE vec_ind = 0;
		for (INTE_TYPE col_ind = p; col_ind < n; col_ind++, X_col += n)
			for (INTE_TYPE row_ind = col_ind + 1; row_ind < n; row_ind++)
				X_col[row_ind] = x[vec_ind++];
		X->mat2vec();
		X->vec2mat();
		dexpSkewSymm_inverse(*Y, *X, *dexp);

		vec_ind = 0;
		for (INTE_TYPE col_ind = p; col_ind < n; col_ind++, Y_col += n)
			for (INTE_TYPE row_ind = col_ind + 1; row_ind < n; row_ind++)
				r[vec_ind++] -= Y_col[row_ind];
	}
};

class Stiefel_Point_Newton
{
	typedef StiefelGeoNewtonOperator Operator;

public:
	INTE_TYPE n;
	INTE_TYPE m;
	INTE_TYPE p;
	INTE_TYPE dm;

	SpecOrthMat Q;
	View_ColMat<REAL_TYPE> Qp;
	View_ColMat<REAL_TYPE> Qm;

	SkewSymmMat S;
	View_ColMat<REAL_TYPE> A;
	View_ColMat<REAL_TYPE> B;
	View_ColMat<REAL_TYPE> C;

	Operator OP;
	SkewSymmMat X;
	SkewSymmMat Y;
	LowerTraversal col_tra_n;
	LowerTraversal blk_tra_n;
	dexpSkewSymmPara dexp_S;

	SkewSymmMat Z;
	LowerTraversal col_tra_m;
	ColMat<REAL_TYPE> **Q_Z;

	matExpPadeApproximant EPA_Z;
	REAL_TYPE *Work_SAF;
	REAL_TYPE *Work_MatMul;

	Stiefel_Point_Newton() : n(0), p(0), m(0), dm(0), Q(SpecOrthMat()), Qp(View_ColMat<REAL_TYPE>()), Qm(View_ColMat<REAL_TYPE>()),
							 S(SkewSymmMat()), A(View_ColMat<REAL_TYPE>()), B(View_ColMat<REAL_TYPE>()), C(View_ColMat<REAL_TYPE>()),
							 OP(Operator()), X(SkewSymmMat()), Y(SkewSymmMat()), col_tra_n(LowerTraversal()), blk_tra_n(LowerTraversal()),
							 dexp_S(dexpSkewSymmPara()), Z(SkewSymmMat()), col_tra_m(LowerTraversal()), Q_Z(nullptr), EPA_Z(matExpPadeApproximant()),
							 Work_SAF(nullptr), Work_MatMul(nullptr){};
	Stiefel_Point_Newton(INTE_TYPE r, INTE_TYPE c) : n(r), m(r - c), p(c), dm((r - c) * (r - c - 1) / 2),
													 Q(SpecOrthMat(r)), S(SkewSymmMat(r)), OP(Operator()), X(SkewSymmMat(r)), Y(SkewSymmMat(r)),
													 col_tra_n(LowerTraversal(r, true, strict_lower_col_traversal)),
													 blk_tra_n(LowerTraversal(r, true, strict_lower_blk_traversal)),
													 col_tra_m(LowerTraversal(r - c, true, strict_lower_col_traversal)),
													 dexp_S(dexpSkewSymmPara(r)), Z(SkewSymmMat(r - c)),
													 Q_Z(new ColMat<REAL_TYPE> *[_StGeoNewton_ExpZ_Scale_Order + 1]), EPA_Z(matExpPadeApproximant(r - c)),
													 Work_SAF(new REAL_TYPE[(r + 2) * r]), Work_MatMul(new REAL_TYPE[r * r])
	{
		Q.fast_col_access();
		Qp = View_ColMat<REAL_TYPE>(Q.v, r, r, c);
		Qm = View_ColMat<REAL_TYPE>(Q.fcol(c), r, r, r - c);

		S.initial_saf(Work_SAF);
		S.fast_col_access();
		S.set_low_vec(&col_tra_n);
		A = View_ColMat<REAL_TYPE>(S.v, r, c, c);
		B = View_ColMat<REAL_TYPE>(S.v + c, r, r - c, c);
		C = View_ColMat<REAL_TYPE>(S.fcol(c) + c, r, r - c, r - c);
		dexp_S.setupSAF(&S.saf);

		OP.n = n;
		OP.p = p;
		OP.dexp = &dexp_S;
		OP.X = &X;
		OP.Y = &Y;
		X.set_low_vec(&col_tra_n);
		Y.set_low_vec(&blk_tra_n);
		Z.set_low_vec(&col_tra_m);

		for (INTE_TYPE i = 0; i < _StGeoNewton_ExpZ_Scale_Order + 1; i++)
			Q_Z[i] = new ColMat<REAL_TYPE>(m, m);
	};
	void copy(const Stiefel_Point_Newton &src)
	{
		this->Q.copy(src.Q);
		this->S.copy(src.S);
		this->dexp_S.copy(src.dexp_S);
		// Shallow objects (View_ColMat) do not participate the deep copy.
		// Workspace X, Y, Z, Q_Z are not copied.
		// Workspace SCF_R is not copied
		// Workspace EPA_Z is not copied
		// Workspace W is not copied.
	};
	Stiefel_Point_Newton(const Stiefel_Point_Newton &src) : Stiefel_Point_Newton(src.n, src.p) { this->copy(src); };
	void swap(Stiefel_Point_Newton &rhs)
	{
		Q.swap(rhs.Q);
		Qp.swap(rhs.Qp);
		Qm.swap(rhs.Qm);
		S.swap(rhs.S);
		A.swap(rhs.A);
		B.swap(rhs.B);
		C.swap(rhs.C);

		this->X.swap(rhs.X);
		this->Y.swap(rhs.Y);
		this->col_tra_n.swap(rhs.col_tra_n);
		this->blk_tra_n.swap(rhs.blk_tra_n);
		this->dexp_S.swap(rhs.dexp_S);
		OP.X = &X;
		OP.Y = &Y;
		OP.dexp = &dexp_S;

		this->Z.swap(rhs.Z);
		this->col_tra_m.swap(rhs.col_tra_m);

		using std::swap;
		swap(this->n, rhs.n);
		swap(this->p, rhs.p);
		swap(this->m, rhs.m);
		swap(this->dm, rhs.dm);
		swap(this->Work_SAF, rhs.Work_SAF);
		swap(this->Work_MatMul, rhs.Work_MatMul);
		swap(this->Q_Z, rhs.Q_Z);
	}
	Stiefel_Point_Newton &operator=(const Stiefel_Point_Newton &rhs)
	{
		Stiefel_Point_Newton temp(rhs);
		swap(temp);
		return (*this);
	}
	~Stiefel_Point_Newton()
	{
		delete[] Work_SAF;
		Work_SAF = nullptr;
		delete[] Work_MatMul;
		Work_MatMul = nullptr;
		for (INTE_TYPE i = 0; i < _StGeoNewton_ExpZ_Scale_Order + 1; i++)
			delete Q_Z[i];
		delete[] Q_Z;
		Q_Z = nullptr;
	};

	void init(const ColMat<REAL_TYPE> &Q0)
	{
		Q.assign(Q0);
		memcpy(Work_MatMul, Q0.v, sizeof(REAL_TYPE) * n * p);
		memset(X.v, 0, sizeof(REAL_TYPE) * n * n);
		X.mat2vec();
		memset(Y.v, 0, sizeof(REAL_TYPE) * n * n);
		Y.mat2vec();
	}

	void Q2Saf(const REAL_TYPE *v)
	{
		S.saf.SchurAngular_SpecOrth(v, n);
		// S.saf.printf();
	};
	void Saf2S()
	{
		S.saf.compute_SkewSymmMat(S);
	};
	void Saf2Dexp()
	{
		dexp_S.setupPara();
		// dexp_S.printf(3);
	}
	void Z2QZ()
	{
		EPA_Z.expPadeApprox(Z, _StGeoNewton_ExpZ_Scale_Order);
		auto ptr = Q_Z;
		(*ptr)->assign(*EPA_Z.Ms[11]);
		for (INTE_TYPE i = 0; i < _StGeoNewton_ExpZ_Scale_Order; i++, ptr++)
			my_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, 1.0, (*ptr)->v, m, (*ptr)->v, m, 0.0, (*(ptr + 1))->v, m);
	};
	void QZ2M(INTE_TYPE i)
	{
		my_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, m, m, 1.0, Qm.v, n, (*(Q_Z + i))->v, m, 0.0, Work_MatMul + n * p, n);
	};
	void M2Q()
	{
		memcpy(Q.v + n * p, Work_MatMul + n * p, sizeof(REAL_TYPE) * n * m);
	};
	INTE_TYPE S2Z(INTE_TYPE &MaxIter, INTE_TYPE Restart, REAL_TYPE &RelTol, REAL_TYPE MaxNorm, ColMat<REAL_TYPE> &H, ColMat<REAL_TYPE> &V, REAL_TYPE *Work)
	{
		INTE_TYPE info;
		S.mat2vec();
		memset(Z.lv, 0, sizeof(REAL_TYPE) * dm);
		info = GMRES(OP, Z.lv, S.lv + col_tra_n.m2v[p + 1 + p * n], H, dm, MaxIter, Restart, RelTol, MaxNorm, Work, V);
		// OP.X->printf("MatX in GMRES:\n");
		// OP.Y->printf("MatY in GMRES:\n");
		my_dscal(dm, -1, Z.lv);
		Z.vec2mat();
		return info;
	};
	REAL_TYPE objval()
	{
		return normsqFro(C);
	};
};