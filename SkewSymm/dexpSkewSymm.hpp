#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

#include "schAngFac.hpp"
#include "skewSymmMat.hpp"
#include "dgemm.hpp"
#include "dgemv.hpp"

class dexpSkewSymmPara
{
    typedef SchurAngularFactor FACTOR;

public:
    // Native members owned by dexpSkewSymmPara:
    REAL_TYPE *forward; // m (m - 1) / 2 * 16 + m * 4 parameters of dexp stored in strict lower 4 x 4 block-column-major order
    REAL_TYPE *inverse; // m (m - 1) / 2 * 16 + m * 4 parameters of dlog stored in strict lower 4 x 4 block-column-major order
    REAL_TYPE *work;    // d x d workspace, for the matrix congruence, i.e., the dgemm for the matrix multiplication.
    // Above members need to be deep-copied in the assignment operator and released in the destructor.

    // Friend member not owne by dexpSkewSymmPara:
    FACTOR *saf;
    REAL_TYPE *a;
    REAL_TYPE *v;
    // Above members only require shallow copy in the assignment operator and they are released by destructor(s) in other class(es).
    // The integrated class SchurAngularFactor is recommended but not necessary.
    // The Schur vectors v and the angles a stored in a continuous memory suffices the task.

    INTE_TYPE d;
    INTE_TYPE asize;
    INTE_TYPE _22_blk_size;
    INTE_TYPE _parvec_size;
    INTE_TYPE _matvec_size;
    INTE_TYPE _lowvec_size;

    dexpSkewSymmPara()
    {
        forward = nullptr;
        inverse = nullptr;
        d = 0;
        asize = 0;
        _22_blk_size = 0;
        _parvec_size = 0;
        _matvec_size = 0;
        _lowvec_size = 0;
    };
    dexpSkewSymmPara(INTE_TYPE dim) : d(dim), asize(dim / 2),
                                      saf(nullptr), a(nullptr), v(nullptr)
    {
        _parvec_size = asize * (asize - 1) * 8 + asize * 4;
        _22_blk_size = (asize * (asize - 1)) / 2;
        _lowvec_size = (d * (d - 1)) / 2;
        _matvec_size = d * d;
        forward = new REAL_TYPE[_parvec_size];
        inverse = new REAL_TYPE[_parvec_size];
        work = new REAL_TYPE[_matvec_size];
    };
    void copy(const dexpSkewSymmPara &src)
    {
        memcpy(this->forward, src.forward, sizeof(REAL_TYPE) * _parvec_size);
        memcpy(this->inverse, src.inverse, sizeof(REAL_TYPE) * _parvec_size);
        // memcpy(this->work, src.work, sizeof(REAL_TYPE) * work);       // Workspace is not copied.
        this->saf = src.saf;
        this->a = this->saf->a;
        this->v = this->saf->v;
    };
    dexpSkewSymmPara(const dexpSkewSymmPara &src) : dexpSkewSymmPara(src.d)
    {
        copy(src);
    };
    dexpSkewSymmPara(FACTOR *factor) : dexpSkewSymmPara(factor->d)
    {
        saf = factor;
        v = factor->v;
        a = factor->a;
    };
    dexpSkewSymmPara(INTE_TYPE dim, REAL_TYPE *schur_vector, REAL_TYPE *angles) : dexpSkewSymmPara(dim)
    {
        saf = nullptr;
        v = schur_vector;
        a = angles;
    };
    void swap(dexpSkewSymmPara &src)
    {
        using std::swap;
        swap(this->forward, src.forward);
        swap(this->inverse, src.inverse);
        swap(this->work, src.work);
        swap(this->saf, src.saf);
        swap(this->a, src.a);
        swap(this->v, src.v);

        swap(this->d, src.d);
        swap(this->asize, src.asize);
        swap(this->_22_blk_size, src._22_blk_size);
        swap(this->_parvec_size, src._parvec_size);
        swap(this->_matvec_size, src._matvec_size);
        swap(this->_lowvec_size, src._lowvec_size);
    };
    dexpSkewSymmPara &operator=(const dexpSkewSymmPara &src)
    {
        dexpSkewSymmPara temp(src);
        swap(temp);
        return (*this);
    };

    ~dexpSkewSymmPara()
    {
        if (forward)
            delete[] forward;
        forward = nullptr;
        if (inverse)
            delete[] inverse;
        inverse = nullptr;
        saf = nullptr;
        v = nullptr;
        a = nullptr;
    }

    void setupWork(INTE_TYPE dim)
    {
        if (d != dim)
        {
            dexpSkewSymmPara temp(dim);
            swap(temp);
        }
    }

    void setupSAF(FACTOR *factor)
    {
        saf = factor;
        v = saf->v;
        a = saf->a;
    }

    void setupSAF(INTE_TYPE dim, REAL_TYPE *schur_vectors, REAL_TYPE *angles)
    {
        setupWork(dim);
        v = schur_vectors;
        a = angles;
    }

    REAL_TYPE _dexpSkewSymm_sxdx(REAL_TYPE x) { return ((x < 1e-15) && (x > -1e-15)) ? (1.0 - x * x / 6.0) : sin(x) / x; };
    REAL_TYPE _dexpSkewSymm_cxdx(REAL_TYPE x) { return ((x < 1e-15) && (x > -1e-15)) ? (-x / 2.0) : (cos(x) - 1.0) / x; };
    REAL_TYPE _dexpSkewSymm_xctx(REAL_TYPE x) { return ((x < 1e-15) && (x > -1e-15)) ? 3.0 / (3.0 + x * x) : x / tan(x); };

    void _setup_22_paras(REAL_TYPE *_forward, REAL_TYPE *_inverse, REAL_TYPE x, REAL_TYPE y);
    void _setup_12_paras(REAL_TYPE *_forward, REAL_TYPE *_inverse, REAL_TYPE x);

    void setupPara();

    void _dexpSkewSymm_congruence(REAL_TYPE *Y, INTE_TYPE ldy, REAL_TYPE *X, INTE_TYPE ldx, BOOL_TYPE forward_action) const;
    void _dexpSkewSymm_congruence(SkewSymmMat &MatY, SkewSymmMat &MatX, BOOL_TYPE forward_action) const
    {
        _dexpSkewSymm_congruence(MatY.v, d, MatX.v, d, forward_action);
        MatY.setdiag();
        MatY.mat2vec();
    }
    void _dexpSkewSymm_congruence(const View_ColMat<REAL_TYPE> &ViewY, const View_ColMat<REAL_TYPE> &ViewX, BOOL_TYPE forward_action) const
    {
        _dexpSkewSymm_congruence(ViewY.v, ViewY.ld, ViewX.v, ViewX.ld, forward_action);
    };

    void _dexpSkewSymm_forward_core(REAL_TYPE *lvY, REAL_TYPE *lvX) const;
    void _dexpSkewSymm_inverse_core(REAL_TYPE *lvY, REAL_TYPE *lvX) const;
    void _dexpSkewSymm_forward_core(SkewSymmMat &Y, SkewSymmMat &X) const;
    void _dexpSkewSymm_inverse_core(SkewSymmMat &Y, SkewSymmMat &X) const;

    void printf(INTE_TYPE subsystem_num)
    {
        INTE_TYPE print_cnt = subsystem_num < _22_blk_size ? subsystem_num : _22_blk_size;
        View_ColMat<REAL_TYPE> sys = View_ColMat<REAL_TYPE>();
        using std::printf;
        printf("\nLinear operator L_{Theta} and its inverse determined by the angles Theta of skew symmetric matrix S.,\n");
        printf("such that Dexp_S[X] = exp(S) * R * (L_{Theta}(R^T * X * R)) * R^{T}.\n");
        printf("The operators consists of independent 4 x 4 systems L_{i,j} acting on i,j-th 2 x 2 matrix partitions\n");
        printf("and the optitional independent 2 x 2 systems L_{i,j} acting on i,j-th 1 x 2 matrix partitions.\n");
        auto forward_ptr = forward;
        auto inverse_ptr = inverse;
        for (INTE_TYPE blk_ind = 0; blk_ind < print_cnt; blk_ind++, forward_ptr += 16, inverse_ptr += 16)
        {
            sys.reset_array(forward_ptr, 4, 4, 4);
            printf("Forward system:\n");
            sys.printf();
            sys.reset_array(inverse_ptr, 4, 4, 4);
            printf("Inverse system:\n");
            sys.printf();
        }
        if (print_cnt < _22_blk_size)
        {
            printf("... other 4 x 4 systems omitted ...\n");
        }
        if (d % 2)
        {
            print_cnt = subsystem_num < asize ? subsystem_num : asize;
            forward_ptr = forward + 16 * _22_blk_size;
            inverse_ptr = inverse + 16 * _22_blk_size;
            for (INTE_TYPE blk_ind = 0; blk_ind < print_cnt; blk_ind++, forward_ptr += 4, inverse_ptr += 4)
            {
                sys.reset_array(forward_ptr, 2, 2, 2);
                printf("Forward system:\n");
                sys.printf();
                sys.reset_array(inverse_ptr, 2, 2, 2);
                printf("Inverse system:\n");
                sys.printf();
            }
            if (print_cnt < asize)
            {
                printf("... other 2 x 2 systems omitted ...\n");
            }
        }
    };
    void printf() { printf(_22_blk_size); };

    friend void dexpSkewSymm_forward(SkewSymmMat &Y, SkewSymmMat &X, dexpSkewSymmPara &Para);
    friend void dexpSkewSymm_inverse(SkewSymmMat &Y, SkewSymmMat &X, dexpSkewSymmPara &Para);
};

// void DexpSkewSymm(SkewSymm::SkewSymmMat &Y, const SkewSymm::DexpPara &para, const SkewSymm::SkewSymmMat &X, REAL_TYPE *work, BOOL_TYPE forward);

// void DexpSkewSymm(SkewSymm::SkewSymmMat &Y, const SkewSymm::DexpPara &para, const REAL_TYPE *vectors, const SkewSymm::SkewSymmMat &X, BOOL_TYPE forward);
