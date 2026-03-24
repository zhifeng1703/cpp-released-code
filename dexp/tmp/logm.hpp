#pragma once

#include "arrVec.hpp"
#include "colMat.hpp"
#include "mmlt.hpp"
#include "matOp.hpp"
#include "skewSchFac.hpp"

class NearbyLogm
{
	// This class is the basic unit of the nearby matrix logarithm on special orthogonal matrices.
	// It keeps (1) the special orthogonal matrix (2) the skew-symmetric matrix (3) the ``shared'' Schur decomposition.
	// For a given reference A, it can perform the following tasks.
	// 		(1) Compute the trial solution X aiming for \|A-X\|_2 < \pi.
	//		(2) Compute the actual distance with the trial solution d = \|A-X\|_2
	//		(3) Verify if A and X are connected.
	// 		(4) In the case where \nlog_A(Q) is not defined, keep necessary information of Q for later calls.
	// A copy of SSF of Q is explicitly stored (the R and A), which can be reused with no cost.
	// Necessary External Objects: (1) SkewSchurFactor SSF, (2) n x n Workspace.

	typedef ColMat<REAL_TYPE> MATX_TYPE;
	typedef ColMat<INTE_TYPE> MATI_TYPE;
	typedef ArrVec<REAL_TYPE> VECT_TYPE;
	typedef SkewSchurFactor FACT_TYPE;
	typedef NearbyLogm SELF_TYPE;

public:
	INTE_TYPE d;
	INTE_TYPE m;

	MATX_TYPE Q;
	MATX_TYPE X;
	MATX_TYPE R;

	VECT_TYPE Arctan;
	VECT_TYPE ShiftA;

	MATI_TYPE LabelI;
	MATI_TYPE LabelJ;

	static constexpr REAL_TYPE TWO_PI = 2.0 * M_PI;

	NearbyLogm() : d(0), m(0), Q(MATX_TYPE()), X(MATX_TYPE()), R(MATX_TYPE()), Arctan(VECT_TYPE()), ShiftA(VECT_TYPE()), LabelI(MATI_TYPE()), LabelJ(MATI_TYPE()) {};
	NearbyLogm(INTE_TYPE dim) : d(dim), m(dim / 2),
								Q(MATX_TYPE(dim, dim)), X(MATX_TYPE(dim, dim)), R(MATX_TYPE(dim, dim)),
								Arctan(VECT_TYPE(dim / 2)), ShiftA(VECT_TYPE(dim / 2)),
								LabelI(MATI_TYPE(dim / 2, dim / 2)), LabelJ(MATI_TYPE(dim / 2, dim / 2)) {};
	void copy(const SELF_TYPE &src)
	{
		this->d = src.d;
		this->m = src.m;
		this->Q.copy(src.Q);
		this->X.copy(src.R);
		this->R.copy(src.R);

		this->Arctan.copy(src.Arctan);
		this->ShiftA.copy(src.ShiftA);

		this->LabelI.copy(src.LabelI);
		this->LabelJ.copy(src.LabelJ);
	}
	NearbyLogm(const SELF_TYPE &src) : NearbyLogm(src.d) { this->copy(src); };
	void swap(SELF_TYPE &rhs)
	{
		std::swap(this->d, rhs.d);
		std::swap(this->m, rhs.m);
		this->Q.swap(rhs.Q);
		this->X.swap(rhs.X);
		this->R.swap(rhs.R);

		this->Arctan.swap(rhs.Arctan);
		this->ShiftA.swap(rhs.ShiftA);

		this->LabelI.swap(rhs.LabelI);
		this->LabelJ.swap(rhs.LabelJ);
	}

	NearbyLogm(SELF_TYPE &&rhs) noexcept : d(rhs.d), m(rhs.m),
										   Q(std::move(rhs.Q)), X(std::move(rhs.X)), R(std::move(rhs.R)),
										   Arctan(std::move(rhs.Arctan)), ShiftA(std::move(rhs.ShiftA)),
										   LabelI(std::move(rhs.LabelI)), LabelJ(std::move(rhs.LabelJ)) {};
	NearbyLogm &operator=(SELF_TYPE rhs)
	{
		swap(rhs);
		return (*this);
	}
	~NearbyLogm() {};

	void Factor_SpecOrth(FACT_TYPE &SSF)
	{
		SSF.Factor_SpecOrth(Q.v, d);
		Arctan.Assign(SSF.A);
		R.Assign(SSF.R.v, SSF.d);
	}

	void Shifting_Delta(REAL_TYPE *Delta, INTE_TYPE ldd, REAL_TYPE *MatA, INTE_TYPE lda, REAL_TYPE *Work, INTE_TYPE ldw)
	{
		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, d, d, d, 1.0, R.v, d, MatA, lda, 0.0, Work, ldw);
		skewblas_mmlt(CblasNoTrans, CblasNoTrans, d, d, d, 1.0, Work, ldw, R.v, d, 0.0, Delta, d);
		skewl2m(Delta, ldd, d);

		ShiftA.Assign(Arctan);

		REAL_TYPE angular_diff = 0;
		INTE_TYPE shift = 0;
		for (auto ind = 0; ind < m; ind++)
		{
			angular_diff = Delta[(2 * ind + 1) + (2 * ind) * ldd] - ShiftA[ind];
			shift = std::lround(angular_diff / TWO_PI);
			ShiftA[ind] += shift * TWO_PI;
			Delta[(2 * ind + 1) + (2 * ind) * ldd] -= ShiftA[ind];
			Delta[(2 * ind) + (2 * ind + 1) * ldd] = -Delta[(2 * ind + 1) + (2 * ind) * ldd];
		}
	};

	REAL_TYPE Delta_Norm(REAL_TYPE *Delta, INTE_TYPE ldd, FACT_TYPE &SSF)
	{
		SSF.Factor_SkewSymm(Delta, ldd);
		if (SSF.a == 0)
			return 0;
		auto ind = cblas_idamax(SSF.a, SSF.A, 1);
		return abs(SSF.A[ind]);
	};

	static REAL_TYPE _gap2nearest(REAL_TYPE s)
	{

		REAL_TYPE t = s / TWO_PI; // target in integer units
		// nearest integer to t
		auto n = std::lround(t);

		if (n == 0)
		{
			// compare distance to +1 and -1 in integer units
			REAL_TYPE dpos = std::fabs(1.0 - t);
			REAL_TYPE dneg = std::fabs(-1.0 - t);
			return std::fabs(((dpos < dneg) ? 1.0 : -1.0) * TWO_PI - s);
		}
		else
			return std::fabs((REAL_TYPE)n * TWO_PI - s);
	}

	REAL_TYPE _dist2locus(INTE_TYPE dim, const REAL_TYPE *VecA)
	{
		INTE_TYPE count = dim / 2;

		REAL_TYPE best = 1e10;
		for (auto i = 0; i < count; i++)
		{
			for (auto j = i + 1; j < count; j++)
			{
				REAL_TYPE s1 = VecA[i] + VecA[j];
				REAL_TYPE gap1 = _gap2nearest(s1) * 0.5;
				if (gap1 < best)
					best = gap1;

				REAL_TYPE s2 = VecA[i] - VecA[j];
				REAL_TYPE gap2 = _gap2nearest(s2) * 0.5;
				if (gap2 < best)
					best = gap2;
			}
		}
		if (dim != 2 * count)
		{
			for (auto i = 0; i < count; i++)
			{
				REAL_TYPE gap = _gap2nearest(VecA[i]);
				if (gap < best)
					best = gap;
			}
		}
		return best;
	}

	REAL_TYPE _dist2locus(const REAL_TYPE *VecA) { return _dist2locus(d, VecA); };

	void Angular_Label(INTE_TYPE dim, const REAL_TYPE *VecA, INTE_TYPE *MatI, INTE_TYPE ldi, INTE_TYPE *MatJ, INTE_TYPE ldj)
	{
		INTE_TYPE count = dim / 2;
		for (auto i = 0; i < m; i++)
		{
			for (auto j = 0; j < m; j++)
			{
				MatI[i + j * ldi] = static_cast<int>(std::trunc((VecA[i] + VecA[j]) / TWO_PI));
				MatJ[i + j * ldj] = static_cast<int>(std::trunc((VecA[i] - VecA[j]) / TWO_PI));
			}
		}
	};

	void Angular_Label() { Angular_Label(d, ShiftA.v, LabelI.v, m, LabelJ.v, m); };
};

class NearbyLogmCurve
{
	// This class is the basic unit of the nearby matrix logarithm on special orthogonal matrices.
	// It keeps (1) the special orthogonal matrix (2) the skew-symmetric matrix (3) the ``shared'' Schur decomposition.
	// For a given reference A, it can perform the following tasks.
	// 		(1) Compute the trial solution X aiming for \|A-X\|_2 < \pi.
	//		(2) Compute the actual distance with the trial solution d = \|A-X\|_2
	//		(3) Verify if A and X are connected.
	// 		(4) In the case where \nlog_A(Q) is not defined, keep necessary information of Q for later calls.
	// A copy of SSF of Q is explicitly stored (the R and A), which can be reused with no cost.
	// Necessary External Objects: (1) SkewSchurFactor SSF, (2) n x n Workspace.

	typedef ColMat<REAL_TYPE> MATX_TYPE;
	typedef ColMat<INTE_TYPE> MATI_TYPE;
	typedef ArrVec<REAL_TYPE> VECT_TYPE;
	typedef SkewSchurFactor FACT_TYPE;
	typedef NearbyLogm SELF_TYPE;

public:
	INTE_TYPE d;
	INTE_TYPE m;

	MATX_TYPE Q;
	MATX_TYPE X;
	MATX_TYPE R;

	VECT_TYPE Arctan;
	VECT_TYPE ShiftA;

	MATI_TYPE LabelI;
	MATI_TYPE LabelJ;

	static constexpr REAL_TYPE TWO_PI = 2.0 * M_PI;

	NearbyLogm() : d(0), m(0), Q(MATX_TYPE()), X(MATX_TYPE()), R(MATX_TYPE()), Arctan(VECT_TYPE()), ShiftA(VECT_TYPE()), LabelI(MATI_TYPE()), LabelJ(MATI_TYPE()) {};
	NearbyLogm(INTE_TYPE dim) : d(dim), m(dim / 2),
								Q(MATX_TYPE(dim, dim)), X(MATX_TYPE(dim, dim)), R(MATX_TYPE(dim, dim)),
								Arctan(VECT_TYPE(dim / 2)), ShiftA(VECT_TYPE(dim / 2)),
								LabelI(MATI_TYPE(dim / 2, dim / 2)), LabelJ(MATI_TYPE(dim / 2, dim / 2)) {};
	void copy(const SELF_TYPE &src)
	{
		this->d = src.d;
		this->m = src.m;
		this->Q.copy(src.Q);
		this->X.copy(src.R);
		this->R.copy(src.R);

		this->Arctan.copy(src.Arctan);
		this->ShiftA.copy(src.ShiftA);

		this->LabelI.copy(src.LabelI);
		this->LabelJ.copy(src.LabelJ);
	}
	NearbyLogm(const SELF_TYPE &src) : NearbyLogm(src.d) { this->copy(src); };
	void swap(SELF_TYPE &rhs)
	{
		std::swap(this->d, rhs.d);
		std::swap(this->m, rhs.m);
		this->Q.swap(rhs.Q);
		this->X.swap(rhs.X);
		this->R.swap(rhs.R);

		this->Arctan.swap(rhs.Arctan);
		this->ShiftA.swap(rhs.ShiftA);

		this->LabelI.swap(rhs.LabelI);
		this->LabelJ.swap(rhs.LabelJ);
	}

	NearbyLogm(SELF_TYPE &&rhs) noexcept : d(rhs.d), m(rhs.m),
										   Q(std::move(rhs.Q)), X(std::move(rhs.X)), R(std::move(rhs.R)),
										   Arctan(std::move(rhs.Arctan)), ShiftA(std::move(rhs.ShiftA)),
										   LabelI(std::move(rhs.LabelI)), LabelJ(std::move(rhs.LabelJ)) {};
	NearbyLogm &operator=(SELF_TYPE rhs)
	{
		swap(rhs);
		return (*this);
	}
	~NearbyLogm() {};

	void Factor_SpecOrth(FACT_TYPE &SSF)
	{
		SSF.Factor_SpecOrth(Q.v, d);
		Arctan.Assign(SSF.A);
		R.Assign(SSF.R.v, SSF.d);
	}

	void Shifting_Delta(REAL_TYPE *Delta, INTE_TYPE ldd, REAL_TYPE *MatA, INTE_TYPE lda, REAL_TYPE *Work, INTE_TYPE ldw)
	{
		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, d, d, d, 1.0, R.v, d, MatA, lda, 0.0, Work, ldw);
		skewblas_mmlt(CblasNoTrans, CblasNoTrans, d, d, d, 1.0, Work, ldw, R.v, d, 0.0, Delta, d);
		skewl2m(Delta, ldd, d);

		ShiftA.Assign(Arctan);

		REAL_TYPE angular_diff = 0;
		INTE_TYPE shift = 0;
		for (auto ind = 0; ind < m; ind++)
		{
			angular_diff = Delta[(2 * ind + 1) + (2 * ind) * ldd] - ShiftA[ind];
			shift = std::lround(angular_diff / TWO_PI);
			ShiftA[ind] += shift * TWO_PI;
			Delta[(2 * ind + 1) + (2 * ind) * ldd] -= ShiftA[ind];
			Delta[(2 * ind) + (2 * ind + 1) * ldd] = -Delta[(2 * ind + 1) + (2 * ind) * ldd];
		}
	};

	REAL_TYPE Delta_Norm(REAL_TYPE *Delta, INTE_TYPE ldd, FACT_TYPE &SSF)
	{
		SSF.Factor_SkewSymm(Delta, ldd);
		if (SSF.a == 0)
			return 0;
		auto ind = cblas_idamax(SSF.a, SSF.A, 1);
		return abs(SSF.A[ind]);
	};

	static REAL_TYPE _gap2nearest(REAL_TYPE s)
	{

		REAL_TYPE t = s / TWO_PI; // target in integer units
		// nearest integer to t
		auto n = std::lround(t);

		if (n == 0)
		{
			// compare distance to +1 and -1 in integer units
			REAL_TYPE dpos = std::fabs(1.0 - t);
			REAL_TYPE dneg = std::fabs(-1.0 - t);
			return std::fabs(((dpos < dneg) ? 1.0 : -1.0) * TWO_PI - s);
		}
		else
			return std::fabs((REAL_TYPE)n * TWO_PI - s);
	}

	REAL_TYPE _dist2locus(INTE_TYPE dim, const REAL_TYPE *VecA)
	{
		INTE_TYPE count = dim / 2;

		REAL_TYPE best = 1e10;
		for (auto i = 0; i < count; i++)
		{
			for (auto j = i + 1; j < count; j++)
			{
				REAL_TYPE s1 = VecA[i] + VecA[j];
				REAL_TYPE gap1 = _gap2nearest(s1) * 0.5;
				if (gap1 < best)
					best = gap1;

				REAL_TYPE s2 = VecA[i] - VecA[j];
				REAL_TYPE gap2 = _gap2nearest(s2) * 0.5;
				if (gap2 < best)
					best = gap2;
			}
		}
		if (dim != 2 * count)
		{
			for (auto i = 0; i < count; i++)
			{
				REAL_TYPE gap = _gap2nearest(VecA[i]);
				if (gap < best)
					best = gap;
			}
		}
		return best;
	}

	REAL_TYPE _dist2locus(const REAL_TYPE *VecA) { return _dist2locus(d, VecA); };

	void Angular_Label(INTE_TYPE dim, const REAL_TYPE *VecA, INTE_TYPE *MatI, INTE_TYPE ldi, INTE_TYPE *MatJ, INTE_TYPE ldj)
	{
		INTE_TYPE count = dim / 2;
		for (auto i = 0; i < m; i++)
		{
			for (auto j = 0; j < m; j++)
			{
				MatI[i + j * ldi] = static_cast<int>(std::trunc((VecA[i] + VecA[j]) / TWO_PI));
				MatJ[i + j * ldj] = static_cast<int>(std::trunc((VecA[i] - VecA[j]) / TWO_PI));
			}
		}
	};

	void Angular_Label() { Angular_Label(d, ShiftA.v, LabelI.v, m, LabelJ.v, m); };
};
