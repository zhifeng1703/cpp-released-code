#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include "hhMat.hpp"
#include "skewMat.hpp"
#include "mmlt.hpp"

class SkewSchurFactor
{
    typedef SkewSchurFactor SELF_TYPE;
    typedef ColMat<REAL_TYPE> MATX_TYPE;
    typedef ArrVec<REAL_TYPE> VECT_TYPE;
    typedef HouseholderMatrix HHMT_TYPE;

public:
    INTE_TYPE d; // Dimension d x d
    INTE_TYPE m; // m = floor(d / 2)
    INTE_TYPE k; // k = ceil(d / 2)
    INTE_TYPE a; // # of nonzero angles a \leq d. This quantity has been disabled and now uses m as default.

    MATX_TYPE R;
    VECT_TYPE A;

    INTE_TYPE o;         // Determinant
    ArrVec<INTE_TYPE> E; // Integer shifts to principal angles

    HHMT_TYPE H;
    MATX_TYPE U;
    MATX_TYPE Vt;

    MATX_TYPE Work;

    SkewSchurFactor() : d(0), m(0), k(0), o(0), R(ColMat<REAL_TYPE>()), A(ArrVec<REAL_TYPE>()), E(ArrVec<INTE_TYPE>()), H(HouseholderMatrix()), U(ColMat<REAL_TYPE>()), Vt(ColMat<REAL_TYPE>()), Work(ColMat<REAL_TYPE>()) {};
    SkewSchurFactor(INTE_TYPE dim) : d(dim), m(dim / 2), k(dim - dim / 2), o(0), R(ColMat<REAL_TYPE>(dim, dim)), A(ArrVec<REAL_TYPE>(dim / 2)), E(ArrVec<INTE_TYPE>(dim / 2)),
                                     H(HouseholderMatrix(dim)),
                                     U(ColMat<REAL_TYPE>(dim / 2, dim / 2)), Vt(ColMat<REAL_TYPE>(dim - dim / 2, dim - dim / 2)),
                                     Work(ColMat<REAL_TYPE>(dim, dim)) {};
    void copy(const SELF_TYPE &src)
    {
        this->d = src.d;
        this->m = src.m;
        this->k = src.k;
        this->o = src.o;

        this->R.copy(src.R);
        this->A.copy(src.A);
        this->E.copy(src.E);

        this->H.copy(src.H);
        this->U.copy(src.U);
        this->Vt.copy(src.Vt);
        // Work space is not copied.
    };
    SkewSchurFactor(const SELF_TYPE &src) : SkewSchurFactor(src.d) { this->copy(src); };
    void swap(SELF_TYPE &rhs)
    {
        this->R.swap(rhs.R);
        this->A.swap(rhs.A);
        this->E.swap(rhs.E);

        this->H.swap(rhs.H);
        this->U.swap(rhs.U);
        this->Vt.swap(rhs.Vt);

        this->Work.swap(rhs.Work);
        using std::swap;
        swap(this->d, rhs.d);
        swap(this->m, rhs.m);
        swap(this->k, rhs.k);
        swap(this->o, rhs.o);
    }
    SkewSchurFactor &operator=(const SELF_TYPE &rhs)
    {
        SELF_TYPE temp(rhs);
        swap(temp);
        return (*this);
    }
    ~SkewSchurFactor() {};

    void _chessboard_scatter(INTE_TYPE dim, REAL_TYPE *board, INTE_TYPE ldb, REAL_TYPE *U, INTE_TYPE ldu, REAL_TYPE *Vt, INTE_TYPE ldvt);

    void _givens_elimin_bidiag_extra_column(REAL_TYPE *_cos, REAL_TYPE *_sin, REAL_TYPE *diag, REAL_TYPE *offd, REAL_TYPE extra);

    void _givens_rotate_cols(INTE_TYPE coli, INTE_TYPE colj, REAL_TYPE _cos, REAL_TYPE _sin, REAL_TYPE *MatVt, INTE_TYPE n);

    void _givens_rotate(INTE_TYPE n, REAL_TYPE *coli, INTE_TYPE inci, REAL_TYPE *colj, INTE_TYPE incj, REAL_TYPE _cos, REAL_TYPE _sin)
    {
        cblas_dcopy(n, coli, inci, Work.v, 1);
        cblas_daxpby(n, _sin, colj, incj, _cos, coli, inci);
        cblas_daxpby(n, -_sin, Work.v, 1, _cos, colj, incj);
    }

    void Factor_SkewSymm(REAL_TYPE *MatS, INTE_TYPE lds)
    {
        H.Assign(MatS, lds);
        SchurAngular_SkewSymm();
    };

    void _principal_angles(REAL_TYPE *MatQ, INTE_TYPE ldq, REAL_TYPE *MatR, INTE_TYPE ldr)
    {
        // Based on the computed Schur vectors in H.MatR and the original matrix,
        // determine the principal angles of the orthogonal matrix in O(n^2).
        REAL_TYPE _cos, _sin, _vec;
        REAL_TYPE *_col_r = MatR;
        a = 0;
        for (auto blk_ind = 0; blk_ind < m; blk_ind++, _col_r += 2 * ldr)
            for (auto row_ind = 0; row_ind < d; row_ind++)
                if (abs(_col_r[row_ind]) > SBLAS_ZERO)
                {
                    _sin = A[blk_ind];
                    _vec = cblas_ddot(d, MatQ + row_ind, ldq, _col_r, 1) - _sin * _col_r[row_ind + d];
                    _cos = std::copysign(sqrt(1 - _sin * _sin), _vec / _col_r[row_ind]);
                    // std::printf("\nsine: %1.8e,\t cosine: %1.8e, \t entry 1: %1.8e, \t entry 2: %1.8e\n", _sin, _cos, _vec, _col_r[row_ind]);
                    A[blk_ind] = atan2(_sin, _cos);
                    a += (abs(A[blk_ind]) > SBLAS_ZERO);
                    break;
                }

        a = m; // disable a for accuracy concern in the Stiefel geodesic problem.
    };

    void _principal_angles(REAL_TYPE *MatQ, INTE_TYPE ldq) { _principal_angles(MatQ, ldq, R.v, d); };

    void _ssf_principal_angles(REAL_TYPE *VecA)
    {
        auto Av = A.v;
        auto Ev = E.v;
        for (auto blk_ind = 0; blk_ind < m; blk_ind++)
        {
            int x = static_cast<int>(std::floor((VecA[blk_ind] + M_PI) / (2 * M_PI)));
            Ev[blk_ind] = x;
            Av[blk_ind] = VecA[blk_ind] - 2 * M_PI * x;

            // safety against roundoff drift outside (-pi, pi]
            if (Av[blk_ind] <= -M_PI)
            {
                Av[blk_ind] += 2 * M_PI;
                Ev[blk_ind] -= 1;
            }
            else if (Av[blk_ind] > M_PI)
            {
                Av[blk_ind] -= 2 * M_PI;
                Ev[blk_ind] += 1;
            }
        }
    };

    void Factor_SpecOrth(REAL_TYPE *MatQ, INTE_TYPE ldq)
    {
        H.Skew_Assign(MatQ, ldq);
        SchurAngular_SkewSymm();
        Explict_Vector();
        _principal_angles(MatQ, ldq);
    }

    void Factor(REAL_TYPE *MatM, INTE_TYPE ldm, BOOL_TYPE isSkew = true)
    {
        if (isSkew)
            Factor_SkewSymm(MatM, ldm);
        else
            Factor_SpecOrth(MatM, ldm);
    };

    inline INTE_TYPE _ssf_sgn(REAL_TYPE val) { return (std::abs(val) < SBLAS_ZERO) ? 0 : (0 < val) - (val < 0); };
    INTE_TYPE _ssf_det(REAL_TYPE *Q, INTE_TYPE ldq)
    {
        if (d == 0)
            return 1;

        Work.Assign(Q, ldq);

        ArrVec<INTE_TYPE> ipiv(d);

        INTE_TYPE info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, d, d, Work.v, d, ipiv.v);
        if (info != 0)
            return 0;

        INTE_TYPE det = 1;

        // Each row interchange flips determinant sign.
        // In LAPACK, ipiv[i] != i+1 means row i+1 was swapped.
        for (INTE_TYPE i = 0; i < d; i++)
            if (ipiv[i] != i + 1)
                det = -det;

        // Multiply by the signs of diagonal entries of U.
        for (int i = 0; i < d; i++)
        {
            det *= _ssf_sgn(Work.v[i + i * d]);
        }
        return det;
    };

    void _ssf_swap_cols(INTE_TYPE i, INTE_TYPE j);
    void _ssf_flip_cols(INTE_TYPE j);

    void SchurAngular_SkewSymm();
    void Canonical_SchurAngular();

    void Explict_Vector(REAL_TYPE *MatR, INTE_TYPE ldr);
    void Explict_Vector() { Explict_Vector(R.v, R.r); };

    void SkewCongruence_Explict(CBLAS_TRANSPOSE trans, REAL_TYPE *MatS, INTE_TYPE lds);
    // void SkewCongruence_Compact(CBLAS_TRANSPOSE trans, REAL_TYPE *MatS, INTE_TYPE lds);
    // void SkewCongruence_Raw(CBLAS_TRANSPOSE trans, REAL_TYPE *MatS, INTE_TYPE lds);

    void SkewCongruence(CBLAS_TRANSPOSE trans, REAL_TYPE *MatS, INTE_TYPE lds)
    {
        SkewCongruence_Explict(trans, MatS, lds);
    };

    void GetSkewSymm(REAL_TYPE *MatA, INTE_TYPE lda)
    {
        Work.Zero();
        for (auto blk_ind = 0; blk_ind < m; blk_ind++)
        {
            cblas_daxpy(d, -A[blk_ind], R.v + (2 * blk_ind + 1) * d, 1, Work.v + 2 * blk_ind, d);
            cblas_daxpy(d, A[blk_ind], R.v + 2 * blk_ind * d, 1, Work.v + 2 * blk_ind + 1, d);
        }
        skewblas_mmlt(CblasNoTrans, CblasNoTrans, d, d, d, 1.0, R.v, d, Work.v, d, 0.0, MatA, lda);
    }

    void Exponential(REAL_TYPE *MatQ, INTE_TYPE ldq)
    {
        R.Copyto(MatQ, ldq);
        for (auto blk_ind = 0; blk_ind < a; blk_ind++)
            this->_givens_rotate(d, MatQ + (2 * blk_ind) * ldq, 1, MatQ + (2 * blk_ind + 1) * ldq, 1, cos(A[blk_ind]), sin(A[blk_ind]));
        Work.Assign(MatQ, ldq);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, d, d, d, 1.0, Work.v, d, R.v, d, 0.0, MatQ, ldq);
    };

    void Exponential(ColMat<REAL_TYPE> &MatQ) { Exponential(MatQ.v, MatQ.r); };
    void Exponential(View_ColMat<REAL_TYPE> &ViewQ) { Exponential(ViewQ.v, ViewQ.ld); };

    void Principal_Logarithm(REAL_TYPE *MatX, INTE_TYPE ldx)
    {
        _ssf_principal_angles(A.v);
        GetSkewSymm(MatX, ldx);
    }

    void Diffeomorphic_Logarithm(REAL_TYPE *MatX, INTE_TYPE ldx, REAL_TYPE *MatS, INTE_TYPE lds, SkewSchurFactor &ssfS)
    {
        // ssfS.Canonical_SchurAngular();
        // this->Canonical_SchurAngular();

        // Both the ssfS and this ssf must be canonically aligned before calling Diffeomorphic_Logarithm. (Note that the cost of one aligment is negligible.)

        const REAL_TYPE two_pi = 2 * M_PI;

        BOOL_TYPE special_component = true;
        if (ssfS.E[0] != 0 && ssfS.E[0] != -1)
            special_component = false;
        else
        {
            for (auto blk_ind = 1; blk_ind < m; blk_ind++)
            {
                if (ssfS.E[blk_ind] != 0)
                {
                    special_component = false;
                    break;
                }
            }
        }

        _ssf_principal_angles(A.v);
        if (!special_component)
        {
            memcpy(E.v, ssfS.E.v, sizeof(INTE_TYPE) * m);
            for (auto blk_ind = 0; blk_ind < m; blk_ind++)
                A.v[blk_ind] += E[blk_ind] * two_pi;
        }
        else
        {
            cblas_dgemv(CblasColMajor, CblasNoTrans, d, d, 1.0, MatS, lds, R.v, 1, 0.0, Work.v, 1);
            REAL_TYPE temp = A.v[0] - cblas_ddot(d, R.v + d, 1, Work.v, 1);
            if (abs(temp) > abs(temp - two_pi))
            {
                E.v[0] = -1;
                A.v[0] -= two_pi;
            }
        }
        GetSkewSymm(MatX, ldx);
    }

    void SpecialDiffeo_Logarithm(REAL_TYPE *MatX, INTE_TYPE ldx, REAL_TYPE *MatS, INTE_TYPE lds)
    {

        const REAL_TYPE two_pi = 2 * M_PI;

        _ssf_principal_angles(A.v);

        cblas_dgemv(CblasColMajor, CblasNoTrans, d, d, 1.0, MatS, lds, R.v, 1, 0.0, Work.v, 1);
        REAL_TYPE temp = A.v[0] - cblas_ddot(d, R.v + d, 1, Work.v, 1);
        if (abs(temp) > abs(temp - two_pi))
        {
            E.v[0] = -1;
            A.v[0] -= two_pi;
        }

        GetSkewSymm(MatX, ldx);
    }

    REAL_TYPE _dist2twopi(REAL_TYPE x)
    {
        const REAL_TYPE two_pi = 2.0 * M_PI;
        INTE_TYPE l = lround(x / two_pi);
        if (l != 0)
            return std::abs(x - l * two_pi);
        else if (x > 0)
            return two_pi - x;
        else
            return two_pi + x;
    }

    REAL_TYPE Dist2ConjugateLocus()
    {
        REAL_TYPE dist = M_PI;
        REAL_TYPE temp = 0;
        if (d % 2)
            for (auto ind = 0; ind < m; ind++)
            {
                temp = _dist2twopi(A.v[ind]);
                if (temp < dist)
                    dist = temp;
            }

        for (auto ind = 0; ind < m; ind++)
            for (auto jnd = ind + 1; jnd < m; jnd++)
            {
                temp = _dist2twopi(A.v[ind] + A.v[jnd]) / 2;
                if (temp < dist)
                    dist = temp;
                temp = _dist2twopi(A.v[ind] - A.v[jnd]) / 2;
                if (temp < dist)
                    dist = temp;
            }

        // A.printf("Angles:\n");
        // std::printf("\nRadius of the inscribed ball:\t %f\n", dist);
        return dist;
    };

#ifndef MATLAB_MEX_BUILD

    void printf(const char *s)
    {
        std::printf("%s", s);
        this->R.printf("Schur Vectors R:\n");
        this->A.printf("Off-diagonals in the block diagonal skew symmetric D where RDR' = A:\n");
    };
    void printf() { this->printf(""); };
#endif
};