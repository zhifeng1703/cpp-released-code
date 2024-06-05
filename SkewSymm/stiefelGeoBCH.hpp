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
#include "blas_lapack_interface/sylvester.hpp"

class Stiefel_Point_BCH
{
public:
	INTE_TYPE n;
	INTE_TYPE m;
	INTE_TYPE p;

	SpecOrthMat Q;
	View_ColMat<REAL_TYPE> Qp;
	View_ColMat<REAL_TYPE> Qm;

	SkewSymmMat S;
	View_ColMat<REAL_TYPE> A;
	View_ColMat<REAL_TYPE> B;
	View_ColMat<REAL_TYPE> C;

	ColMat<REAL_TYPE> R;

	SkewSymmMat Z;
	ColMat<REAL_TYPE> Q_Z;

	SchurCanonicalFactor SCF_R;
	matExpPadeApproximant EPA_Z;
	ColMat<REAL_TYPE> W;
	View_ColMat<REAL_TYPE> Wm;
	View_ColMat<REAL_TYPE> Wmm;

	Stiefel_Point_BCH() : n(0), p(0), m(0), Q(SpecOrthMat()), Qp(View_ColMat<REAL_TYPE>()), Qm(View_ColMat<REAL_TYPE>()),
						  S(SkewSymmMat()), A(View_ColMat<REAL_TYPE>()), B(View_ColMat<REAL_TYPE>()), C(View_ColMat<REAL_TYPE>()),
						  R(ColMat<REAL_TYPE>()), Z(SkewSymmMat()), Q_Z(ColMat<REAL_TYPE>()),
						  SCF_R(SchurCanonicalFactor()), EPA_Z(matExpPadeApproximant()),
						  W(ColMat<REAL_TYPE>()), Wm(View_ColMat<REAL_TYPE>()), Wmm(View_ColMat<REAL_TYPE>()){};

	Stiefel_Point_BCH(INTE_TYPE r, INTE_TYPE c) : n(r), m(r - c), p(c), Q(SpecOrthMat(r)), S(SkewSymmMat(r)),
												  R(ColMat<REAL_TYPE>(r - c, r - c)), Q_Z(ColMat<REAL_TYPE>(r - c, r - c)), Z(SkewSymmMat(r - c)),
												  SCF_R(SchurCanonicalFactor(r - c)), EPA_Z(matExpPadeApproximant(r - c)), W(ColMat<REAL_TYPE>(r, r))
	{

		Q.initial_saf(W.v);
		Q.fast_col_access();
		Qp = View_ColMat<REAL_TYPE>(Q.v, r, r, c);
		Qm = View_ColMat<REAL_TYPE>(Q.fcol(c), r, r, r - c);

		S.fast_col_access();
		A = View_ColMat<REAL_TYPE>(S.v, r, c, c);
		B = View_ColMat<REAL_TYPE>(S.v + c, r, r - c, c);
		C = View_ColMat<REAL_TYPE>(S.fcol(c) + c, r, r - c, r - c);

		Z.fast_col_access();

		W.fast_col_access();
		Wm = View_ColMat<REAL_TYPE>(W.fcol(c), r, r, r - c);
		Wmm = View_ColMat<REAL_TYPE>(W.fcol(c), r, r - c, r - c);
	};
	void copy(const Stiefel_Point_BCH &src)
	{
		this->Q.copy(src.Q);
		this->S.copy(src.S);
		this->R.copy(src.R);
		this->Q_Z.copy(src.Q_Z);
		this->Z.copy(src.Z);
		// Shallow objects (View_ColMat) do not participate the deep copy.
		// Workspace SCF_R is not copied
		// Workspace EPA_Z is not copied
		// Workspace W is not copied.
	};
	Stiefel_Point_BCH(const Stiefel_Point_BCH &src) : Stiefel_Point_BCH(src.n, src.p) { this->copy(src); };
	void swap(Stiefel_Point_BCH &rhs)
	{
		Q.swap(rhs.Q);
		Qp.swap(rhs.Qp);
		Qm.swap(rhs.Qm);
		S.swap(rhs.S);
		A.swap(rhs.A);
		B.swap(rhs.B);
		C.swap(rhs.C);
		R.swap(rhs.R);
		SCF_R.swap(rhs.SCF_R);
		Q_Z.swap(rhs.Q_Z);
		Z.swap(rhs.Z);
		EPA_Z.swap(rhs.EPA_Z);
		W.swap(rhs.W);
		Wm.swap(rhs.Wm);
		Wmm.swap(rhs.Wmm);

		using std::swap;
		swap(this->n, rhs.n);
		swap(this->p, rhs.p);
		swap(this->m, rhs.m);
	}
	Stiefel_Point_BCH &operator=(const Stiefel_Point_BCH &rhs)
	{
		Stiefel_Point_BCH temp(rhs);
		swap(temp);
		return (*this);
	}

	~Stiefel_Point_BCH(){};

	void Update()
	{
		S2Z();
		Z2Q();
		Q2S();
	};

	void Q2S()
	{
		Q.SchurAngular_SpecOrth();

		// printf("Schur Factor of Q, vectors:\n");
		// Q.saf.svec.printf();
		// printf("Schur Factor of Q, angles:\n");
		// for (INTE_TYPE a_ind = 0; a_ind < Q.saf.nzsize; a_ind++)
		//	printf("%1.3f\t", Q.saf.a[a_ind]);
		// printf("\n");

		Q.saf.compute_SkewSymmMat(S);
	};
	void Z2QZ()
	{
		EPA_Z.expPadeApprox(Z);
		EPA_Z.exp(Q_Z);
	};
	void QZ2Q()
	{
		my_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, m, m, 1.0, Qm.v, n, Q_Z.v, m, 0.0, Wm.v, n);
		memcpy(Qm.v, Wm.v, sizeof(REAL_TYPE) * n * m);
	};
	void Z2Q()
	{
		Z2QZ();
		QZ2Q();
	};
	void S2Z()
	{
		memset(R.v, 0, sizeof(REAL_TYPE) * m * m);
		for (INTE_TYPE ind = 0; ind < m; ind++)
			R.v[ind + ind * m] = -0.5;
		my_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, m, p, 1.0 / 12.0, B.v, B.ld, B.v, B.ld, 1.0, R.v, m);

		SCF_R.setup_decompose_matrix(R);
		SCF_R.factor();

		symsyl(1, Z, R, C, SCF_R, Wmm);
		// scal(-1, Z);
		Z.Skew();
	};

	REAL_TYPE objval()
	{
		return normsqFro(C);
	};
};