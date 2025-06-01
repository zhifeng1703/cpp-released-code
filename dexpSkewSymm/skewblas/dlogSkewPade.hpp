#pragma once

#include "mkl.h"
#include "colMat.hpp"
#include "skewSchFac.hpp"
#include "pivotedLU.hpp"

// This code implements the scaling and squaring method of the matrix exponential proposed in [Higham09].
// Not that it is the olde version that does not contained special treatment of the scaling order.
// For the improved version of the scaling and squaring mehtod, see expmss_improved.hpp.

// [Higham09]: Higham, Nicholas J. "The scaling and squaring method for the matrix exponential revisited."
// SIAM review 51, no. 4 (2009): 747-764.

// Note that this code serve as the preprocessing step (in the worst case scenario) of the directional
// derivative dexp implemented in dexpPadeSeries, which requires [13/13] Pade approximant.
// Therefore, this code always assumes the [13/13] Pade approximant is in use, without preprocessing on
// reducing the matrix 1-norm, i.e., it implements line 17 - 21 in [Higham09](Algorithm 2.3)

// #define _EXPM_PADE_APPROX_MAX_SCALING 10
// #define _EXPM_PADE_APPROX_MAX_PARANUM 13

// const INTE_TYPE _EXPM_PADE_APPROX_3[4] = {120, 60, 12, 1};
// const INTE_TYPE _EXPM_PADE_APPROX_5[6] = {30240, 15120, 3360, 420, 30, 1};
// const INTE_TYPE _EXPM_PADE_APPROX_7[8] = {17297280, 8652600, 1995840, 277200, 25200, 1512, 56, 1};
// const INTE_TYPE _EXPM_PADE_APPROX_9[10] = {17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1};
// const INTE_TYPE _EXPM_PADE_APPROX_13[14] = {64764752532480000, 32382376266240000, 7771770303897600,
//                                             1187353796428800, 129060195264000, 10559470521600,
//                                             670442572800, 33522128640, 1323241920,
//                                             40840800, 960960, 16380, 182, 1};

// const REAL_TYPE _EXPM_PADE_APPROX_BOUNDS[13] = {2.11e-8, 3.56e-4, 1.08e-2, 6.49e-2, 2.00e-1, 4.37e-1, 7.83e-1, 1.23e0, 1.78e0, 2.42e0, 3.13e0, 3.90e0, 4.74e0};

const REAL_TYPE _LOGM_PADE_ALPHA_7[] = {0.20897959, 0.19091503, 0.19091503, 0.13985270, 0.13985270, 0.06474248, 0.06474248};
const REAL_TYPE _LOGM_PADE_ALPHA_6[] = {0.23395697, 0.23395697, 0.18038079, 0.18038079, 0.08566225, 0.08566225};
const REAL_TYPE _LOGM_PADE_ALPHA_5[] = {0.28444444, 0.23931434, 0.23931434, 0.11846344, 0.11846344};
const REAL_TYPE _LOGM_PADE_ALPHA_4[] = {0.17392742, 0.17392742, 0.32607258, 0.32607258};
const REAL_TYPE _LOGM_PADE_ALPHA_3[] = {0.44444444, 0.27777778, 0.27777778};
const REAL_TYPE _LOGM_PADE_ALPHA_2[] = {0.50000000, 0.50000000};
const REAL_TYPE _LOGM_PADE_ALPHA_1[] = {1.00000000};
const REAL_TYPE *_LOGM_PADE_ALPHA[] = {_LOGM_PADE_ALPHA_1, _LOGM_PADE_ALPHA_2, _LOGM_PADE_ALPHA_3, _LOGM_PADE_ALPHA_4, _LOGM_PADE_ALPHA_5, _LOGM_PADE_ALPHA_6, _LOGM_PADE_ALPHA_7};

const REAL_TYPE _LOGM_PADE_BETA_7[] = {0.50000000, 0.70292258, 0.29707742, 0.87076559, 0.12923481, 0.97455396, 0.02544604};
const REAL_TYPE _LOGM_PADE_BETA_6[] = {0.61930959, 0.38069041, 0.83060469, 0.16939531, 0.96623476, 0.03376524};
const REAL_TYPE _LOGM_PADE_BETA_5[] = {0.50000000, 0.76923466, 0.23076534, 0.95308992, 0.04691008};
const REAL_TYPE _LOGM_PADE_BETA_4[] = {0.93056816, 0.06943184, 0.66999052, 0.33000948};
const REAL_TYPE _LOGM_PADE_BETA_3[] = {0.50000000, 0.88729833, 0.11270167};
const REAL_TYPE _LOGM_PADE_BETA_2[] = {0.78867513, 0.21132487};
const REAL_TYPE _LOGM_PADE_BETA_1[] = {0.50000000};
const REAL_TYPE *_LOGM_PADE_BETA[] = {_LOGM_PADE_BETA_1, _LOGM_PADE_BETA_2, _LOGM_PADE_BETA_3, _LOGM_PADE_BETA_4, _LOGM_PADE_BETA_5, _LOGM_PADE_BETA_6, _LOGM_PADE_BETA_7};

const REAL_TYPE _LOGM_PADE_APPROX_BOUNDS[] = {2.11e-8, 3.56e-4, 1.08e-2, 6.49e-2, 2.00e-1, 4.37e-1, 7.83e-1, 1.23e0, 1.78e0, 2.42e0, 3.13e0, 3.90e0, 4.74e0};

class dlogSkewPadeApprox
{
    typedef dlogSkewPadeApprox SELF_TYPE;
    typedef ArrVec<REAL_TYPE> VECT_TYPE;
    typedef ColMat<REAL_TYPE> MATX_TYPE;

public:
    INTE_TYPE d;
    INTE_TYPE s;
    INTE_TYPE m;

    const REAL_TYPE *alpha;
    const REAL_TYPE *beta;
    // Pre-computed constant type coefficients, not managed by this object.

    VECT_TYPE A;
    MATX_TYPE R;

    MATX_TYPE X;
    MATX_TYPE Y;

    MATX_TYPE Work;

    dlogSkewPadeApprox() : s(0), d(0), m(0), alpha(nullptr), beta(nullptr), A(VECT_TYPE()), R(MATX_TYPE()), X(MATX_TYPE()), Y(MATX_TYPE()), Work(MATX_TYPE()) {};
    dlogSkewPadeApprox(INTE_TYPE dim) : s(0), d(dim), m(0),
                                        A(VECT_TYPE(dim / 2)), R(MATX_TYPE(dim, dim)), X(MATX_TYPE(dim, dim)), Y(MATX_TYPE(dim, dim)), Work(MATX_TYPE(dim, dim + 1)) {};
    void copy(const SELF_TYPE &src)
    {
        // assert(d == src.d);
        s = src.s;
        m = src.m;
        this->alpha = src.alpha;
        this->beta = src.beta;
        this->A.copy(src.A);
        this->R.copy(src.R);
        this->X.copy(src.X);
        this->Y.copy(src.Y);
    }
    dlogSkewPadeApprox(const SELF_TYPE &src) : SELF_TYPE(src.d) { this->copy(src); };
    void swap(SELF_TYPE &src)
    {
        this->A.swap(src.A);
        this->R.swap(src.R);
        this->X.swap(src.X);
        this->Y.swap(src.Y);
        this->Work.swap(src.Work);
        using std::swap;
        swap(this->d, src.d);
        swap(this->s, src.s);
        swap(this->m, src.m);
        swap(this->alpha, src.alpha);
        swap(this->beta, src.beta);
    }
    SELF_TYPE &operator=(const SELF_TYPE &rhs)
    {
        SELF_TYPE temp = SELF_TYPE(rhs);
        swap(temp);
        return (*this);
    }
    ~dlogSkewPadeApprox() {};

    void Assign(REAL_TYPE *MatA, INTE_TYPE lda) { X.Assign(MatA, lda); };
    void Assign(const View_ColMat<REAL_TYPE> &ViewA) { Assign(ViewA.v, ViewA.ld); };
    void Assign(const ColMat<REAL_TYPE> &MatA) { Assign(MatA.v, MatA.r); };

    // void _get_orders()
    // {
    //     REAL_TYPE norm = VecMatM[0]->Norm1();
    //     if (norm < _EXPM_PADE_APPROX_BOUNDS[2])
    //     {
    //         m = 3;
    //         _vmsize = 7;
    //         s = 0;
    //     }
    //     else if (norm < _EXPM_PADE_APPROX_BOUNDS[4])
    //     {
    //         m = 5;
    //         _vmsize = 8;
    //         s = 0;
    //     }
    //     else if (norm < _EXPM_PADE_APPROX_BOUNDS[6])
    //     {
    //         m = 7;
    //         _vmsize = 9;
    //         s = 0;
    //     }
    //     else if (norm < _EXPM_PADE_APPROX_BOUNDS[8])
    //     {
    //         m = 9;
    //         _vmsize = 10;
    //         s = 0;
    //     }
    //     else
    //     {
    //         m = 13;
    //         _vmsize = 13;
    //         s = ceil(log2(norm / _EXPM_PADE_APPROX_BOUNDS[12]));
    //         s = s < 0 ? 0 : s;
    //     }
    //     _vqsize = s + 1;
    // };

    void _set_orders(INTE_TYPE order, INTE_TYPE scale)
    {
        s = scale;
        m = order;
        alpha = _LOGM_PADE_ALPHA[m - 1];
        beta = _LOGM_PADE_BETA[m - 1];
    };

    void _set_R(REAL_TYPE *MatQ, INTE_TYPE ldq) { R.Assign(MatQ, ldq); };
    void _set_A(REAL_TYPE *VecA) { A.Assign(VecA); };

    void Parameter(REAL_TYPE *MatQ, INTE_TYPE ldq, REAL_TYPE *VecA)
    {
        _set_R(MatQ, ldq);
        _set_A(VecA);
    }
    void Parameter(const SkewSchurFactor &SSF) { this->Parameter(SSF.R.v, d, SSF.A.v); };
    void Parameter(ColMat<REAL_TYPE> &MatQ, ArrVec<REAL_TYPE> &VecA) { this->Parameter(MatQ.v, d, VecA.v); };

    // void _solve_sylvester(INTE_TYPE scale_order)
    // {
    //     REAL_TYPE *TmpY = X.v, *TmpX = Y.v, *TmpLeft = Work.v;
    //     REAL_TYPE *TmpRight = (d > 3) ? TmpLeft + 4 : TmpLeft;
    //     INTE_TYPE blk_dim = A.d;
    //     REAL_TYPE _scale = 0, _temp = 0;
    //     for (auto col_ind = 0; col_ind < blk_dim; col_ind++)
    //     {
    //         TmpRight[0] = 1 + cos(A.v[col_ind] / pow(2, scale_order));
    //         TmpRight[1] = sin(A.v[col_ind] / pow(2, scale_order));
    //         TmpRight[2] = -TmpRight[1];
    //         TmpRight[3] = TmpRight[0];

    //         TmpY = X.v + (2 * col_ind) * d;

    //         for (auto row_ind = 0; row_ind < blk_dim; row_ind++)
    //         {
    //             TmpLeft[0] = 1 + cos(A.v[row_ind] / pow(2, scale_order));
    //             TmpLeft[1] = sin(A.v[row_ind] / pow(2, scale_order));
    //             TmpLeft[2] = -TmpLeft[1];
    //             TmpLeft[3] = TmpLeft[0];

    //             TmpX[0] = TmpY[0];
    //             TmpX[1] = TmpY[1];
    //             TmpX[2] = TmpY[d];
    //             TmpX[3] = TmpY[d + 1];

    //             LAPACKE_dtrsyl(LAPACK_COL_MAJOR, 'N', 'N', 1, 2, 2, TmpLeft, 2, TmpRight, 2, TmpX, 2, &_scale);

    //             TmpY[0] = TmpX[0] / _scale;
    //             TmpY[1] = TmpX[1] / _scale;
    //             TmpY[d] = TmpX[2] / _scale;
    //             TmpY[d + 1] = TmpX[3] / _scale;

    //             TmpY += 2;
    //         }
    //     }

    //     if (2 * blk_dim != d)
    //     {
    //         for (auto ind = 0; ind < blk_dim; ind++)
    //         {
    //             _scale = tan(A.v[ind] / pow(2, s + 1)) * 0.5;
    //             _temp = X(2 * ind, d - 1);
    //             X(2 * ind, d - 1) = 0.5 * _temp - _scale * X(2 * ind + 1, d - 1);
    //             X(2 * ind + 1, d - 1) = _scale * _temp + 0.5 * X(2 * ind + 1, d - 1);

    //             _temp = X(d, 2 * ind);
    //             X(d - 1, 2 * ind) = 0.5 * _temp + _scale * X(d - 1, 2 * ind + 1);
    //             X(d - 1, 2 * ind + 1) = -_scale * _temp + 0.5 * X(d - 1, 2 * ind + 1);
    //         }
    //     }
    // }

    void _solve_sylvester(INTE_TYPE scale_order)
    {
        REAL_TYPE *currX = X.v;
        INTE_TYPE blk_dim = A.d;
        REAL_TYPE _cos_i, _sin_i, _cos_j, _sin_j;
        REAL_TYPE _ac, _bc, _ad, _bd;
        REAL_TYPE _scale = 0, _temp = 0;
        for (auto col_ind = 0; col_ind < blk_dim; col_ind++, currX = X.v + (2 * col_ind) * d)
            for (auto row_ind = 0; row_ind < blk_dim; row_ind++, currX += 2)
            {
                // Wrong implementation
                _cos_i = cos(A.v[row_ind]);
                _sin_i = sin(A.v[row_ind]);
                _cos_j = cos(A.v[col_ind]);
                _sin_j = sin(A.v[col_ind]);

                _scale = 0.5 / (_cos_i + _cos_j);
                _ac = _scale * (1.0 + _cos_i * _cos_j);
                _bc = _scale * (_sin_i * _cos_j);
                _ad = _scale * (_cos_i * _sin_j);
                _bd = _scale * (_sin_i * _sin_j);

                Work.v[0] = currX[0];
                Work.v[1] = currX[1];
                Work.v[2] = currX[d];
                Work.v[3] = currX[d + 1];

                currX[0] = _ac * Work.v[0] + _bc * Work.v[1] - _ad * Work.v[2] - _bd * Work.v[3];
                currX[1] = -_bc * Work.v[0] + _ac * Work.v[1] + _bd * Work.v[2] - _ad * Work.v[3];
                currX[d] = _ad * Work.v[0] + _bd * Work.v[1] + _ac * Work.v[2] + _bc * Work.v[3];
                currX[d + 1] = -_bd * Work.v[0] + _ad * Work.v[1] - _bc * Work.v[2] + _ac * Work.v[3];
            }

        if (2 * blk_dim != d)
        {
            for (auto ind = 0; ind < blk_dim; ind++)
            {
                _scale = tan(A.v[ind] / pow(2, s + 1)) * 0.5;
                _temp = X(2 * ind, d - 1);
                X(2 * ind, d - 1) = 0.5 * _temp - _scale * X(2 * ind + 1, d - 1);
                X(2 * ind + 1, d - 1) = _scale * _temp + 0.5 * X(2 * ind + 1, d - 1);

                _temp = X(d, 2 * ind);
                X(d - 1, 2 * ind) = 0.5 * _temp + _scale * X(d - 1, 2 * ind + 1);
                X(d - 1, 2 * ind + 1) = -_scale * _temp + 0.5 * X(d - 1, 2 * ind + 1);
            }
        }
    }

    // void _solve_sylvester(INTE_TYPE scale_order)
    // {

    //     INTE_TYPE blk_dim = A.d;
    //     REAL_TYPE _scale = 0, _temp = 0;

    //     Y.Zero();
    //     REAL_TYPE *diagY = Y.v;

    //     for (auto ind = 0; ind < blk_dim; ind++, diagY += 2 + 2 * d)
    //     {
    //         diagY[0] = 1 + cos(A.v[ind] / pow(2, scale_order));
    //         diagY[1] = sin(A.v[ind] / pow(2, scale_order));
    //         diagY[d] = -diagY[1];
    //         diagY[d + 1] = diagY[0];
    //     }
    //     if (d != 2 * blk_dim)
    //         diagY[0] = 1;

    //     LAPACKE_dtrsyl(LAPACK_COL_MAJOR, 'N', 'N', 1, d, d, Y.v, d, Y.v, d, X.v, d, &_scale);
    //     cblas_dscal(d * d, 1.0 / _scale, X.v, 1);
    // }

    void _solve_system(REAL_TYPE _beta, REAL_TYPE *WorkL, REAL_TYPE *WorkV)
    {
        REAL_TYPE _cos, _a, _b;
        memcpy(WorkL, X.v, sizeof(REAL_TYPE) * d * d);
        for (auto ind = 0; ind < A.d; ind++)
        {
            _cos = cos(A.v[ind] / pow(2, s));
            _b = 2 * (_beta * _beta - 1) * (1 - _cos) + 1;
            _a = (1 - _beta * (1 - _cos)) / _b;
            _b = -_beta * sin(A.v[ind] / pow(2, s)) / _b;

            cblas_dcopy(d, WorkL + 2 * ind, d, WorkV, 1);
            cblas_daxpby(d, -_b, WorkL + 2 * ind + 1, d, _a, WorkL + 2 * ind, d);
            cblas_daxpby(d, _b, WorkV, 1, _a, WorkL + 2 * ind + 1, d);

            cblas_dcopy(d, WorkL + (2 * ind) * d, 1, WorkV, 1);
            cblas_daxpby(d, _b, WorkL + (2 * ind + 1) * d, 1, _a, WorkL + (2 * ind) * d, 1);
            cblas_daxpby(d, -_b, WorkL + (2 * ind) * d, 1, _a, WorkL + (2 * ind + 1) * d, 1);
        }
    }

    void Dlog(REAL_TYPE *MatM, INTE_TYPE ldm)
    {
        Assign(MatM, ldm);
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, d, d, d, 1.0, R.v, d, X.v, d, 0.0, Work.v, d);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, 1.0, Work.v, d, R.v, d, 0.0, X.v, d);

        for (auto ind = 0; ind < s; ind++)
            _solve_sylvester(ind);

        Y.Zero();
        for (auto ind = 0; ind < m; ind++)
        {
            _solve_system(beta[ind], Work.v, Work.v + d * d);
            cblas_daxpy(d * d, alpha[ind], Work.v, 1, Y.v, 1);
        }

        cblas_dscal(d * d, pow(2, s), Y.v, 1);

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, 1.0, R.v, d, Y.v, d, 0.0, Work.v, d);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, d, d, d, 1.0, Work.v, d, R.v, d, 0.0, Y.v, d);
    }

    void Dlog(REAL_TYPE *MatN, INTE_TYPE ldn, REAL_TYPE *MatM, INTE_TYPE ldm)
    {
        Dlog(MatM, ldm);
        Y.Copyto(MatN, ldn);
    }
};