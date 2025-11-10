#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

#include "arrVec.hpp"
#include "skewMat.hpp"
#include "mmlt.hpp"
#include "skewSchFac.hpp"

class dexpSkewSymmPara
{

    typedef dexpSkewSymmPara SELF_TYPE;
    typedef LowerTraversal TRAV_TYPE;
    typedef ArrVec<REAL_TYPE> VECT_TYPE;
    typedef ColMat<REAL_TYPE> MATX_TYPE;

public:
    INTE_TYPE d;
    INTE_TYPE a;
    BOOL_TYPE rank_deficient;
    INTE_TYPE _22_blk_size;
    INTE_TYPE _parvec_size;
    INTE_TYPE _matvec_size;
    INTE_TYPE _lowvec_size;

    VECT_TYPE _forward_para; // m (m - 1) / 2 * 16 + m * 4 = _paravec_size parameters of dexp stored in strict lower 4 x 4 block-column-major order
    VECT_TYPE _inverse_para; // m (m - 1) / 2 * 16 + m * 4 = _paravec_size parameters of dlog stored in strict lower 4 x 4 block-column-major order
    VECT_TYPE A;
    MATX_TYPE R;

    // Above members need to be deep-copied in the assignment operator and released in the destructor.

    TRAV_TYPE _blk_tra;
    MATX_TYPE Work; // d x d workspace, for the matrix congruence, i.e., the dgemm for the matrix multiplication.

    dexpSkewSymmPara() {};
    dexpSkewSymmPara(INTE_TYPE dim)
    {
        d = dim;
        a = dim / 2;
        rank_deficient = true;
        _parvec_size = a * (a - 1) * 8 + a * 4;
        _22_blk_size = (a * (a - 1)) / 2;
        _lowvec_size = (d * (d - 1)) / 2;
        _matvec_size = d * d;

        _forward_para = ArrVec<REAL_TYPE>(_parvec_size);
        _inverse_para = ArrVec<REAL_TYPE>(_parvec_size);
        A = ArrVec<REAL_TYPE>(a);
        R = ColMat<REAL_TYPE>(d, d);

        _blk_tra = LowerTraversal(dim, BLK_OFFD);
        Work = ColMat<REAL_TYPE>(d, d);
    };
    void copy(const SELF_TYPE &src)
    {
        this->rank_deficient = src.rank_deficient;
        this->_forward_para.copy(src._forward_para);
        this->_inverse_para.copy(src._inverse_para);
        this->A.copy(src.A);
        this->R.copy(src.R);
    };
    void swap(SELF_TYPE &src)
    {
        this->_forward_para.swap(src._forward_para);
        this->_inverse_para.swap(src._inverse_para);
        this->A.swap(src.A);
        this->R.swap(src.R);

        this->_blk_tra.swap(src._blk_tra);
        this->Work.swap(src.Work);

        using std::swap;
        swap(this->d, src.d);
        swap(this->a, src.a);
        swap(this->rank_deficient, src.rank_deficient);
        swap(this->_22_blk_size, src._22_blk_size);
        swap(this->_parvec_size, src._parvec_size);
        swap(this->_matvec_size, src._matvec_size);
        swap(this->_lowvec_size, src._lowvec_size);
    }
    dexpSkewSymmPara(const SELF_TYPE &src) : dexpSkewSymmPara(src.d) { this->copy(src); };
    dexpSkewSymmPara &operator=(const SELF_TYPE &src)
    {
        SELF_TYPE temp(src);
        swap(temp);
        return (*this);
    };
    ~dexpSkewSymmPara() {};

    void _set_R(REAL_TYPE *MatQ, INTE_TYPE ldq) { R.Assign(MatQ, ldq); };
    void _set_A(REAL_TYPE *VecA);

    void Parameter(REAL_TYPE *MatQ, INTE_TYPE ldq, REAL_TYPE *VecA)
    {
        _set_R(MatQ, ldq);
        _set_A(VecA);
    }
    void Parameter(const SkewSchurFactor &SSF) { this->Parameter(SSF.R.v, d, SSF.A.v); };
    void Parameter(ColMat<REAL_TYPE> &MatQ, ArrVec<REAL_TYPE> &VecA) { this->Parameter(MatQ.v, d, VecA.v); };

    REAL_TYPE _dexpSkewSymm_sxdx(REAL_TYPE x) { return ((x < 1e-15) && (x > -1e-15)) ? (1.0 - x * x / 6.0) : sin(x) / x; };
    REAL_TYPE _dexpSkewSymm_cxdx(REAL_TYPE x) { return ((x < 1e-15) && (x > -1e-15)) ? (-x * 0.5) : (cos(x) - 1.0) / x; };
    REAL_TYPE _dexpSkewSymm_xctx(REAL_TYPE x) { return ((x < 1e-15) && (x > -1e-15)) ? 3.0 / (3.0 - x * x) : x / tan(x); };

    void _setup_22_paras(REAL_TYPE *_forward, REAL_TYPE *_inverse, REAL_TYPE x, REAL_TYPE y);
    void _setup_12_paras(REAL_TYPE *_forward, REAL_TYPE *_inverse, REAL_TYPE x);

    void _dexpSkewSymm_forward_core(REAL_TYPE *lvY, REAL_TYPE *lvX);
    void _dexpSkewSymm_inverse_core(REAL_TYPE *lvY, REAL_TYPE *lvX);
    void _dexpSkewSymm_forward_core(SkewSymmMat &Y, SkewSymmMat &X);
    void _dexpSkewSymm_inverse_core(SkewSymmMat &Y, SkewSymmMat &X);

    void printf(const char *s, INTE_TYPE subsystem_num)
    {
        std::printf("%s", s);
        INTE_TYPE print_cnt = (subsystem_num < _22_blk_size) ? subsystem_num : _22_blk_size;
        View_ColMat<REAL_TYPE> sys = View_ColMat<REAL_TYPE>();
        using std::printf;
        this->R.printf("Schur Vectors:\n");
        this->A.printf("Angles:\n");
        printf("\nLinear operator L_{Theta} and its inverse determined by the angles Theta of skew symmetric matrix S.,\n");
        printf("such that Dexp_S[X] = exp(S) * R * (L_{Theta}(R^T * X * R)) * R^{T}.\n");
        printf("The operators consists of independent 4 x 4 systems L_{i,j} acting on i,j-th 2 x 2 matrix partitions\n");
        printf("and the optitional independent 2 x 2 systems L_{i,j} acting on i,j-th 1 x 2 matrix partitions.\n");
        auto forward_ptr = _forward_para.v;
        auto inverse_ptr = _inverse_para.v;
        printf("Displaying the first %lld systems in the size of 2 x 2:\n----------------------------------------------------------------\n", print_cnt);

        INTE_TYPE angle_row = 1;
        INTE_TYPE angle_col = 0;
        for (INTE_TYPE blk_ind = 0; blk_ind < print_cnt; blk_ind++, forward_ptr += 16, inverse_ptr += 16)
        {
            sys.Reset(forward_ptr, 4, 4, 4);
            printf("Systems with angle[%lld]= %1.3f,\t and angle[%lld]= %1.3f:\n", angle_row, A[angle_row], angle_col, A[angle_col]);
            printf("Multiples to PI, sum:\t %1.3f, diff:\t %1.3f\n", (A[angle_row] + A[angle_col]) / M_PI, (A[angle_row] - A[angle_col]) / M_PI);

            printf("Forward system\n");
            sys.printf();
            sys.Reset(inverse_ptr, 4, 4, 4);
            printf("Inverse system:\n");
            sys.printf();
            angle_row++;
            if (angle_row == a)
            {
                angle_col++;
                angle_row = angle_col + 1;
            }
        }
        if (print_cnt < _22_blk_size)
            printf("... other 4 x 4 systems omitted ...\n----------------------------------------------------------------\n");
        else
            printf("----------------------------------------------------------------\n");

        if (d % 2)
        {
            print_cnt = subsystem_num < a ? subsystem_num : a;
            forward_ptr = _forward_para.v + 16 * _22_blk_size;
            inverse_ptr = _inverse_para.v + 16 * _22_blk_size;
            printf("Displaying the first %lld 2 x 2 systems:\n----------------------------------------------------------------\n", print_cnt);

            for (INTE_TYPE blk_ind = 0; blk_ind < print_cnt; blk_ind++, forward_ptr += 4, inverse_ptr += 4)
            {
                sys.Reset(forward_ptr, 2, 2, 2);
                printf("Forward system:\n");
                sys.printf();
                sys.Reset(inverse_ptr, 2, 2, 2);
                printf("Inverse system:\n");
                sys.printf();
            }
            if (print_cnt < a)
                printf("... other 2 x 2 systems omitted ...\n----------------------------------------------------------------\n");
            else
                printf("----------------------------------------------------------------\n");
        }
    };

    void printf(const char *s) { printf(s, 2); };
    void printf() { printf("", 2); };

    void SkewCongruence(CBLAS_TRANSPOSE trans, REAL_TYPE *MatY, INTE_TYPE ldy, REAL_TYPE *MatX, INTE_TYPE ldx)
    {
        if (trans == CblasNoTrans)
        {
            // Computing Y <- RXR' in the lower triangular part.
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, d, d, d, 1.0, MatX, ldx, R.v, d, 0.0, Work.v, d);
            // View_ColMat<REAL_TYPE>(Work.v, d, d, d).printf();
            skewblas_mmlt(CblasNoTrans, CblasNoTrans, d, d, d, 1.0, R.v, d, Work.v, d, 0.0, MatY, ldy);
            // View_ColMat<REAL_TYPE>(MatY, ldy, d, d).printf();
        }
        else
        {
            // Computing Y <- R'XR in the lower triangular part.
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, d, d, d, 1.0, R.v, d, MatX, ldx, 0.0, Work.v, d);
            // View_ColMat<REAL_TYPE>(Work.v, d, d, d).printf();
            skewblas_mmlt(CblasNoTrans, CblasNoTrans, d, d, d, 1.0, Work.v, d, R.v, d, 0.0, MatY, ldy);
            // View_ColMat<REAL_TYPE>(MatY, ldy, d, d).printf();
        }
    }

    void Forward(SkewSymmMat &Y, SkewSymmMat &X)
    {
        SkewCongruence(CblasTrans, Y.v, d, X.v, d);
        _dexpSkewSymm_forward_core(Y, Y);
        SkewCongruence(CblasNoTrans, Y.v, d, Y.v, d);
        Y.low2upp();
    }

    void Inverse(SkewSymmMat &Y, SkewSymmMat &X)
    {
        SkewCongruence(CblasTrans, Y.v, d, X.v, d);
        _dexpSkewSymm_inverse_core(Y, Y);
        SkewCongruence(CblasNoTrans, Y.v, d, Y.v, d);
        Y.low2upp();
    };
};
