#pragma once

#include "colMat.hpp"
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

#define _EXPM_PADE_APPROX_MAX_SCALING 10
#define _EXPM_PADE_APPROX_MAX_PARANUM 13

const INTE_TYPE _EXPM_PADE_APPROX_3[4] = {120, 60, 12, 1};
const INTE_TYPE _EXPM_PADE_APPROX_5[6] = {30240, 15120, 3360, 420, 30, 1};
const INTE_TYPE _EXPM_PADE_APPROX_7[8] = {17297280, 8652600, 1995840, 277200, 25200, 1512, 56, 1};
const INTE_TYPE _EXPM_PADE_APPROX_9[10] = {17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1};
const INTE_TYPE _EXPM_PADE_APPROX_13[14] = {64764752532480000, 32382376266240000, 7771770303897600,
                                            1187353796428800, 129060195264000, 10559470521600,
                                            670442572800, 33522128640, 1323241920,
                                            40840800, 960960, 16380, 182, 1};

const REAL_TYPE _EXPM_PADE_APPROX_BOUNDS[13] = {2.11e-8, 3.56e-4, 1.08e-2, 6.49e-2, 2.00e-1, 4.37e-1, 7.83e-1, 1.23e0, 1.78e0, 2.42e0, 3.13e0, 3.90e0, 4.74e0};

class expmPadeApprox
{
    typedef expmPadeApprox SELF_TYPE;
    typedef ColMat<REAL_TYPE> MATX_TYPE;
    typedef ArrVec<ColMat<REAL_TYPE> *> VECT_TYPE;
    typedef PivotedLUFactor FACT_TYPE;

public:
    INTE_TYPE d;
    INTE_TYPE s;
    INTE_TYPE m;

    VECT_TYPE VecMatM;
    VECT_TYPE VecMatQ;
    FACT_TYPE PLU;

    INTE_TYPE _vmsize;
    INTE_TYPE _vqsize;

    // [Higham09]: Higham, Nicholas J. "The scaling and squaring method for the matrix exponential revisited."
    // SIAM review 51, no. 4 (2009): 747-764.

    // Order m = 13 (Largest supported order)

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

    expmPadeApprox() : s(0), d(0), m(0), VecMatM(VECT_TYPE()), VecMatQ(VECT_TYPE()), PLU(FACT_TYPE()) {};
    expmPadeApprox(INTE_TYPE dim) : s(0), d(dim), m(0),
                                    VecMatM(VECT_TYPE(_EXPM_PADE_APPROX_MAX_PARANUM)),
                                    VecMatQ(VECT_TYPE(_EXPM_PADE_APPROX_MAX_SCALING + 1)),
                                    PLU(FACT_TYPE(dim)),
                                    _vmsize(0), _vqsize(0) {};
    void copy(const SELF_TYPE &src)
    {
        // assert(d == src.d);
        s = src.s;
        m = src.m;
        _check_and_initialize(this->VecMatM, src._vmsize, src.d);
        _vmsize = src._vmsize;
        for (auto ind = 0; ind < _vmsize; ind++)
            VecMatM[ind]->copy(*src.VecMatM[ind]);
        _check_and_initialize(this->VecMatQ, src._vqsize, src.d);
        _vqsize = src._vqsize;
        for (auto ind = 0; ind < _vqsize; ind++)
            VecMatQ[ind]->copy(*src.VecMatQ[ind]);
        PLU.copy(src.PLU);
    }
    expmPadeApprox(const SELF_TYPE &src) : expmPadeApprox(src.d) { this->copy(src); };
    void swap(SELF_TYPE &src)
    {
        this->VecMatM.swap(src.VecMatM);
        this->VecMatQ.swap(src.VecMatQ);
        this->PLU.swap(src.PLU);
        using std::swap;
        swap(this->d, src.d);
        swap(this->s, src.s);
        swap(this->m, src.m);
        swap(this->_vmsize, src._vmsize);
        swap(this->_vqsize, src._vqsize);
    }
    SELF_TYPE &operator=(const SELF_TYPE &rhs)
    {
        SELF_TYPE temp = SELF_TYPE(rhs);
        swap(temp);
        return (*this);
    }
    ~expmPadeApprox()
    {
        _release_vec_mat(VecMatM);
        _release_vec_mat(VecMatQ);
    };

    void _check_and_initialize(VECT_TYPE &vec, INTE_TYPE length, INTE_TYPE dim)
    {
        for (auto ind = 0; ind < length; ind++)
            if (!vec.v[ind])
                vec.v[ind] = new MATX_TYPE(dim, dim);
            else if (vec.v[ind]->r != dim)
            {
                delete vec.v[ind];
                vec.v[ind] = new MATX_TYPE(dim, dim);
            }
    };

    void _release_vec_mat(VECT_TYPE &vec)
    {
        for (auto ind = 0; ind < vec.d; ind++)
            if (vec.v[ind])
                delete vec.v[ind];
    };

    void Assign(REAL_TYPE *MatA, INTE_TYPE lda) { VecMatM[0]->Assign(MatA, lda); };
    void Assign(const View_ColMat<REAL_TYPE> &ViewA) { Assign(ViewA.v, ViewA.ld); };
    void Assign(const ColMat<REAL_TYPE> &MatA) { Assign(MatA.v, MatA.r); };

    void _get_orders()
    {
        REAL_TYPE norm = VecMatM[0]->Norm1();
        if (norm < _EXPM_PADE_APPROX_BOUNDS[2])
        {
            m = 3;
            _vmsize = 7;
            s = 0;
        }
        else if (norm < _EXPM_PADE_APPROX_BOUNDS[4])
        {
            m = 5;
            _vmsize = 8;
            s = 0;
        }
        else if (norm < _EXPM_PADE_APPROX_BOUNDS[6])
        {
            m = 7;
            _vmsize = 9;
            s = 0;
        }
        else if (norm < _EXPM_PADE_APPROX_BOUNDS[8])
        {
            m = 9;
            _vmsize = 10;
            s = 0;
        }
        else
        {
            m = 13;
            _vmsize = 13;
            s = ceil(log2(norm / _EXPM_PADE_APPROX_BOUNDS[12]));
            s = s < 0 ? 0 : s;
        }
        _vqsize = s + 1;
    };

    void _set_orders(INTE_TYPE order)
    {
        switch (order)
        {
        case 3:
            m = 3;
            _vmsize = 7;
            s = 0;
            break;
        case 5:
            m = 5;
            _vmsize = 8;
            s = 0;
            break;
        case 7:
            m = 7;
            _vmsize = 9;
            s = 0;
            break;
        case 9:
            m = 7;
            _vmsize = 9;
            s = 0;
            break;
        default:
            m = 13;
            _vmsize = 13;
            s = ceil(log2(VecMatM[0]->Norm1() / _EXPM_PADE_APPROX_BOUNDS[12]));
            s = s < 0 ? 0 : s;
            break;
        }
        _vqsize = s + 1;
    };
    void _scale()
    {
        if (s != 0)
            VecMatM[0]->Scale(1.0 / pow(2, s));
    };

    void _square()
    {
        _check_and_initialize(VecMatQ, s + 1, d);
        VecMatQ[0]->Assign(*VecMatM[_vmsize - 1]);
        // VecMatQ[0]->printf("Initial Q_s = exp(A/2^s):\n");
        for (auto ind = 0; ind < s; ind++)
        {
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, 1.0, VecMatQ[ind]->v, d, VecMatQ[ind]->v, d, 0.0, VecMatQ[ind + 1]->v, d);
            // std::printf("Iterating Q_j = exp(A/2^j), j=s-1:-1:0 with j = %i:\n", ind + 1);
            // VecMatQ[ind + 1]->printf();
        }
    };

    void _expm_pade_low();
    void _expm_pade_13();

    void _expm_pade()
    {
        if (m < 13)
            _expm_pade_low();
        else
            _expm_pade_13();
    };

    void Pade(REAL_TYPE *MatA, INTE_TYPE lda)
    {
        _check_and_initialize(VecMatM, 1, d);
        Assign(MatA, lda);
        _get_orders();
        _scale();
    }

    void Expm(REAL_TYPE *MatA, INTE_TYPE lda)
    {
        _check_and_initialize(VecMatM, 1, d);
        Assign(MatA, lda);
        _get_orders();
        _scale();
        _expm_pade();
        _square();
        // VecMatQ[s] stores the exponenant.
    }

    void Expm(REAL_TYPE *MatA, INTE_TYPE lda, INTE_TYPE order)
    {
        _check_and_initialize(VecMatM, 1, d);
        Assign(MatA, lda);
        _set_orders(order);
        _scale();
        _expm_pade();
        _square();
        // VecMatQ[s] stores the exponenant.
    }
    void Expm(MATX_TYPE &MatA) { Expm(MatA.v, MatA.r); };
    void Expm(View_ColMat<REAL_TYPE> &ViewA) { Expm(ViewA.v, ViewA.ld); };
    void Expm(MATX_TYPE &MatA, INTE_TYPE order) { Expm(MatA.v, MatA.r, order); };
    void Expm(View_ColMat<REAL_TYPE> &ViewA, INTE_TYPE order) { Expm(ViewA.v, ViewA.ld, order); };

    void GetExpm(REAL_TYPE *MatQ, INTE_TYPE ldq) { VecMatQ[s]->Copyto(MatQ, ldq); }

    void Expm(REAL_TYPE *MatQ, INTE_TYPE ldq, REAL_TYPE *MatA, INTE_TYPE lda)
    {
        Expm(MatA, lda);
        GetExpm(MatQ, ldq);
    }
    void Expm(MATX_TYPE &MatQ, MATX_TYPE &MatS) { Expm(MatQ.v, MatQ.r, MatS.v, MatS.r); };
    void Expm(View_ColMat<REAL_TYPE> &ViewQ, View_ColMat<REAL_TYPE> &ViewS) { Expm(ViewQ.v, ViewQ.ld, ViewS.v, ViewS.ld); };

    void Expm(REAL_TYPE *MatQ, INTE_TYPE ldq, REAL_TYPE *MatA, INTE_TYPE lda, INTE_TYPE order)
    {
        Expm(MatA, lda, order);
        GetExpm(MatQ, ldq);
    }
    void Expm(MATX_TYPE &MatQ, MATX_TYPE &MatS, INTE_TYPE order) { Expm(MatQ.v, MatQ.r, MatS.v, MatS.r, order); };
    void Expm(View_ColMat<REAL_TYPE> &ViewQ, View_ColMat<REAL_TYPE> &ViewS, INTE_TYPE order) { Expm(ViewQ.v, ViewQ.ld, ViewS.v, ViewS.ld, order); };
};