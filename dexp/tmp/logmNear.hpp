#pragma once

#include "arrVec.hpp"
#include "colMat.hpp"
#include "mmlt.hpp"
#include "matOp.hpp"
#include "skewSchFac.hpp"

class LogmSkewSymm
{
	typedef ColMat<REAL_TYPE> MATX_TYPE;
	typedef ArrVec<REAL_TYPE> VECT_TYPE;
	typedef SkewSchurFactor FACT_TYPE;
	typedef LogmSkewSymm SELF_TYPE;

public:
	INTE_TYPE d;
	INTE_TYPE m;
	MATX_TYPE A;
	MATX_TYPE R;
	MATX_TYPE Delta;
	MATX_TYPE Work;
	VECT_TYPE Theta;
	FACT_TYPE SSF;

	LogmSkewSymm() : d(0), m(0), A(MATX_TYPE()), R(MATX_TYPE()), Delta(MATX_TYPE()), Work(MATX_TYPE()), Theta(VECT_TYPE()), SSF(FACT_TYPE()) {};
	LogmSkewSymm(INTE_TYPE dim) : d(dim), m(dim / 2),
								  A(MATX_TYPE(dim, dim)), R(MATX_TYPE(dim, dim)), Delta(MATX_TYPE(dim, dim)), Work(MATX_TYPE(dim, dim)),
								  Theta(VECT_TYPE(dim / 2)), SSF(FACT_TYPE(dim)) {};
	void copy(const SELF_TYPE &src)
	{
		this->d = src.d;
		this->m = src.m;
		this->A.copy(src.A);
		this->R.copy(src.R);
		this->Delta.copy(src.Delta);
		this->Work.copy(src.Work);
		this->Theta.copy(src.Theta);
		this->SSF.copy(src.SSF);
	}
	LogmSkewSymm(const SELF_TYPE &src) : LogmSkewSymm(src.d) { this->copy(src); };
	void swap(SELF_TYPE &rhs)
	{
		std::swap(this->d, rhs.d);
		std::swap(this->m, rhs.m);
		this->A.swap(rhs.A);
		this->R.swap(rhs.R);
		this->Delta.swap(rhs.Delta);
		this->Work.swap(rhs.Work);
		this->Theta.swap(rhs.Theta);
		this->SSF.swap(rhs.SSF);
	}
	LogmSkewSymm &operator=(const SELF_TYPE &rhs)
	{
		SELF_TYPE temp(rhs);
		swap(temp);
		return (*this);
	}
	~LogmSkewSymm() {};

	void _diff(REAL_TYPE *MatA, INTE_TYPE lda, REAL_TYPE *MatR, INTE_TYPE ldr, REAL_TYPE *Sigma)
	{
		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, d, d, d, 1.0, MatR, ldr, MatA, lda, 0.0, Work.v, d);
		skewblas_mmlt(CblasNoTrans, CblasNoTrans, d, d, d, 1.0, Work.v, d, MatR, ldr, 0.0, Delta.v, d);
		skewl2m(Delta);
		REAL_TYPE angular_diff = 0;
		const REAL_TYPE two_pi = 2.0 * M_PI;
		INTE_TYPE shift = 0;
		for (auto ind = 0; ind < m; ind++)
		{
			angular_diff = Delta(2 * ind + 1, 2 * ind) - Sigma[ind];
			shift = std::lround(angular_diff / two_pi);
			Delta(2 * ind + 1, 2 * ind) -= (Sigma[ind] + shift * two_pi);
			Delta(2 * ind, 2 * ind + 1) = -Delta(2 * ind + 1, 2 * ind);
		}
	};
	REAL_TYPE _norm(SkewSchurFactor &external_ssf)
	{
		external_ssf.Factor_SkewSymm(Delta.v, d);
		REAL_TYPE norm = 0.0;
		for (auto ind = 0; ind < external_ssf.a; ind++)
			if (abs(external_ssf.A[ind]) > norm)
				norm = abs(external_ssf.A[ind]);
		printf("Distance to the reference\t %f:\n", norm);
		return norm;
	}

	BOOL_TYPE _logm_recursive_iteration(REAL_TYPE *MatA, INTE_TYPE lda, REAL_TYPE *MatR, INTE_TYPE ldr, REAL_TYPE *Sigma, REAL_TYPE bound);
};

// #pragma once

// #include "mkl.h"
// #include "dynArr.hpp"
// #include "colMat.hpp"
// #include "dexpSkew.hpp"

// const INTE_TYPE _LOGM_SKEWSYMM_MAXITER = 1000;
// const REAL_TYPE _LOGM_SKEWSYMM_EABSTOL = 1e-12;

// struct _LogmNearRecord
// {
// 	typedef DynamicArr<REAL_TYPE> RECX_TYPE;
// 	typedef DynamicArr<ColMat<REAL_TYPE> *> RECM_TYPE;
// 	typedef DynamicArr<ArrVec<REAL_TYPE> *> RECV_TYPE;

// 	RECX_TYPE Error;
// 	RECX_TYPE NormV;
// 	RECX_TYPE GridT;

// 	RECM_TYPE Q;
// 	RECM_TYPE S;
// 	RECM_TYPE X;
// 	RECM_TYPE Y;
// 	RECM_TYPE R;

// 	RECV_TYPE A;

// 	_LogmNearRecord(INTE_TYPE d, INTE_TYPE c)
// 	{
// 		Error = RECX_TYPE(d, c);
// 		NormV = RECX_TYPE(d, c);
// 		GridT = RECX_TYPE(d, c);

// 		Q = RECM_TYPE(d, c);
// 		S = RECM_TYPE(d, c);
// 		X = RECM_TYPE(d, c);
// 		Y = RECM_TYPE(d, c);
// 		R = RECM_TYPE(d, c);

// 		A = RECV_TYPE(d, c);
// 	};
// 	_LogmNearRecord(INTE_TYPE n) : _LogmNearRecord(n < 50 ? n : 50, (n / 10) < 50 ? 50 : (n / 10)) {};
// 	_LogmNearRecord() : _LogmNearRecord(0, 1) {};

// 	void swap(_LogmNearRecord &src)
// 	{
// 		Error.swap(src.Error);
// 		NormV.swap(src.NormV);
// 		GridT.swap(src.GridT);
// 		Q.swap(src.Q);
// 		S.swap(src.S);
// 		X.swap(src.X);
// 		Y.swap(src.Y);
// 		R.swap(src.R);
// 		A.swap(src.A);
// 	}

// 	~_LogmNearRecord()
// 	{
// 		Q.manual_release_pointers();
// 		S.manual_release_pointers();
// 		X.manual_release_pointers();
// 		Y.manual_release_pointers();
// 		R.manual_release_pointers();
// 		A.manual_release_pointers();
// 	}
// };

// class LogmSkewSymm
// {
// 	typedef LogmSkewSymm SELF_TYPE;
// 	typedef ColMat<REAL_TYPE> MATX_TYPE;
// 	typedef SkewSchurFactor FACT_TYPE;
// 	typedef dexpSkewSymmPara DEXP_TYPE;
// 	typedef _LogmNearRecord RECD_TYPE;

// public:
// 	INTE_TYPE d;
// 	INTE_TYPE iter;

// 	INTE_TYPE _max_iter;

// 	FACT_TYPE SSF;
// 	DEXP_TYPE Dexp;
// 	MATX_TYPE Work;

// 	RECD_TYPE _rec;
// 	BOOL_TYPE _record;
// 	INTE_TYPE _record_container;
// 	INTE_TYPE _record_capacity;

// 	LogmSkewSymm() : d(0), n(0), iter(0), MaxIter(0), ssf(FACT_TYPE()), dexp(DEXP_TYPE()), rec(RECD_TYPE(0, 0)){};
// 	LogmSkewSymm(INTE_TYPE dim, INTE_TYPE maxiter, BOOL_TYPE record) : d(dim), n(0), iter(0), MaxIter(maxiter), ssf(FACT_TYPE(dim)), dexp(DEXP_TYPE(dim))
// 	{
// 		if (record)
// 			rec = RECD_TYPE(maxiter);
// 		else
// 			rec = RECD_TYPE(4, 4);
// 	};
// 	void copy(const SELF_TYPE &src)
// 	{
// 		d = src.d;
// 		n = 0;
// 		MaxIter = src.MaxIter;

// 		ssf = FACT_TYPE(d);
// 		dexp = DEXP_TYPE(d);
// 	}
// 	LogmSkewSymm(const SELF_TYPE &src) : LogmSkewSymm(src.d) { this->copy(src); };
// 	void swap(SELF_TYPE & src)
// 	{
// 		this->ssf.swap(src.ssf);
// 		using std::swap;
// 		swap(this->d, src.d);
// 		swap(this->n, src.n);
// 	}

// 	SELF_TYPE &operator=(const SELF_TYPE &rhs)
// 	{
// 		SELF_TYPE temp = SELF_TYPE(rhs);
// 		swap(temp);
// 		return (*this);
// 	}
// 	~LogmSkewSymm(){};

// 	void _release_records() {
// 	};
// };