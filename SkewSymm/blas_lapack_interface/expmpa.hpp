#pragma once

#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"
#include "matOp.hpp"

#include "dgees.hpp"
#include "pivotedLUFac.hpp"

// This code implements the scaling and squaring method of the matrix exponential proposed in [Higham09].
// Not that it is the olde version that does not contained special treatment of the scaling order.
// For the improved version of the scaling and squaring mehtod, see expmss_improved.hpp.

// [Higham09]: Higham, Nicholas J. "The scaling and squaring method for the matrix exponential revisited."
// SIAM review 51, no. 4 (2009): 747-764.

// Note that this code serve as the preprocessing step (in the worst case scenario) of the directional
// derivative dexp implemented in dexpPadeSeries, which requires [13/13] Pade approximant.
// Therefore, this code always assumes the [13/13] Pade approximant is in use, without preprocessing on
// reducing the matrix 1-norm, i.e., it implements line 17 - 21 in [Higham09](Algorithm 2.3)

#define _EXPM_PADE_MAXIMUM_SCALING_ORDER 16
#define _EXPM_PADE_APPROXIMANT_PARANUM 13

const INTE_TYPE _EXPM_PADE_APPROXIMANT[14] = {64764752532480000, 32382376266240000, 7771770303897600,
                                              1187353796428800, 129060195264000, 10559470521600,
                                              670442572800, 33522128640, 1323241920,
                                              40840800, 960960, 16380, 182, 1};
// const INTE_TYPE _EXPM_PADE_APPROXIMANT[14] = {1, 182, 16380, 960960, 40840800, 1323241920, 33522128640, 670442572800,
//                                                 10559470521600, 129060195264000, 1187353796428800, 7771770303897600,
//                                                 32382376266240000, 64764752532480000};
const REAL_TYPE _EXPM_THETA_COEFFICIENT[14] = {0.0, 3.7e-8, 5.3e-4, 1.5e-2, 8.5e-2, 2.5e-1, 5.4e-1, 9.5e-1, 1.5e0, 2.1e0, 2.8e0, 3.6e0, 4.5e0, 5.4e0};

//{1.495585217958292e-2, 2.539398330063230e-1, 9.504178996162932e-1,
//                                                    2.097847961257068e0, 5.371920351148152e0};

class matExpPadeApproximant
{
public:
    ColMat<REAL_TYPE> **Ms;
    PivotedLUFactor plu;
    INTE_TYPE d;
    INTE_TYPE s;

    // [Higham09]: Higham, Nicholas J. "The scaling and squaring method for the matrix exponential revisited."
    // SIAM review 51, no. 4 (2009): 747-764.

    // Ms = new ColMat<REAL_TYPE> *[_EXPM_PADE_APPROXIMANT_PARANUM];
    // ind:     0       1       2       3       4       5       6       7       8       9       10      11      12
    // Ms       A       A2      A4      A6      W1      W2      Z1      Z2      W       U       V       R       Q

    // A    = X / 2 ^ s is the scaled version
    // A2   = A^2   = A * A
    // A4   = A^4   = A2 * A2
    // A6   = A^6   = A4 * A2
    // W1   = b13 A6 + b11 A4 + b9 A2
    // W2   = b7 A6 + b5 A4 + b3 A2
    // Z1   = b12 A6 + b10 A4 + b8 A2
    // Z2   = b6 A6 + b4 A4 + b2 A2
    // W    = A6 W1 + W2
    // U    = AW
    // V    = A6Z1 + Z2
    // Q    = exp(A)

    matExpPadeApproximant() : Ms(nullptr), plu(PivotedLUFactor()), s(0), d(0){};
    matExpPadeApproximant(INTE_TYPE dim) : Ms(new ColMat<REAL_TYPE> *[_EXPM_PADE_APPROXIMANT_PARANUM]), plu(PivotedLUFactor(dim)), s(0), d(dim)
    {
        for (INTE_TYPE i = 0; i < _EXPM_PADE_APPROXIMANT_PARANUM; i++)
            Ms[i] = new ColMat<REAL_TYPE>(dim, dim);
    };
    matExpPadeApproximant(const matExpPadeApproximant &src) : matExpPadeApproximant(src.d)
    {
        this->s = src.s;
        for (INTE_TYPE i = 0; i < _EXPM_PADE_APPROXIMANT_PARANUM; i++)
            *Ms[i] = *src.Ms[i];
        this->plu = src.plu;
    };
    void swap(matExpPadeApproximant &src)
    {
        this->plu.swap(src.plu);
        using std::swap;
        swap(this->Ms, src.Ms);
        swap(this->d, src.d);
        swap(this->s, src.s);
    }
    matExpPadeApproximant &operator=(const matExpPadeApproximant &rhs)
    {
        matExpPadeApproximant temp = matExpPadeApproximant(rhs);
        swap(temp);
        return (*this);
    }
    ~matExpPadeApproximant()
    {
        if (Ms)
        {
            for (INTE_TYPE i = 0; i < _EXPM_PADE_APPROXIMANT_PARANUM; i++)
                delete Ms[i];
            delete[] Ms;
        }
        Ms = nullptr;
    };

    void copy(const matExpPadeApproximant &src)
    {
        this->s = src.s;
        for (INTE_TYPE i = 0; i < _EXPM_PADE_APPROXIMANT_PARANUM; i++)
            Ms[i]->copy(*src.Ms[i]);
        this->plu.copy(src.plu);
    }

    void setup_perturb(REAL_TYPE *A, INTE_TYPE lda)
    {
        REAL_TYPE *src_ptr = A;
        REAL_TYPE *des_ptr = Ms[0]->v;
        for (INTE_TYPE col_ind = 0; col_ind < d; col_ind++, des_ptr += d, src_ptr += lda)
            memcpy(des_ptr, src_ptr, sizeof(REAL_TYPE) * d);
    };
    void setup_perturb(const View_ColMat<REAL_TYPE> &ViewA)
    {
        setup_perturb(ViewA.v, ViewA.ld);
    };
    void setup_perturb(const ColMat<REAL_TYPE> &A)
    {
        memcpy(Ms[0]->v, A.v, sizeof(REAL_TYPE) * d * d);
    };

    INTE_TYPE scaling_order()
    {
        INTE_TYPE scale = ceil(log2(Ms[0]->opnorm(1) / _EXPM_THETA_COEFFICIENT[13]));
        if (scale < 0)
            return 0;
        else
            return scale;
    };

    void expPadeApprox(INTE_TYPE scaling_order);
    void expPadeApprox() { expPadeApprox(scaling_order()); };
    void expPadeApprox(const ColMat<REAL_TYPE> &A)
    {
        setup_perturb(A);
        expPadeApprox();
    }
    void expPadeApprox(const ColMat<REAL_TYPE> &A, INTE_TYPE scaling_order)
    {
        setup_perturb(A);
        expPadeApprox(scaling_order);
    }
    void expPadeApprox(REAL_TYPE *A, INTE_TYPE lda)
    {
        setup_perturb(A, lda);
        expPadeApprox();
    }
    void exp(ColMat<REAL_TYPE> &A)
    {
        using std::swap;
        // ColMat<REAL_TYPE> *curr_R = &A;
        // ColMat<REAL_TYPE> *next_R = Ms[12];

        ColMat<REAL_TYPE> *curr_R = Ms[12];
        ColMat<REAL_TYPE> *next_R = &A;
        if (s % 2)
        {
            my_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, 1.0, Ms[11]->v, d, Ms[11]->v, d, 0.0, next_R->v, d);
            swap(curr_R, next_R);
            for (INTE_TYPE s_ind = 1; s_ind < s; s_ind++, swap(curr_R, next_R))
                my_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, 1.0, curr_R->v, d, curr_R->v, d, 0.0, next_R->v, d);
        }
        else
        {
            next_R->copy(*Ms[11]);
            swap(curr_R, next_R);
            for (INTE_TYPE s_ind = 0; s_ind < s; s_ind++, swap(curr_R, next_R))
                my_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, 1.0, curr_R->v, d, curr_R->v, d, 0.0, next_R->v, d);
        }
    }
};