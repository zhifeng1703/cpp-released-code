#pragma once

#include <cassert>
#include <cstring>
#include <cmath>

#ifdef MATLAB_MEX_BUILD
#include "blasType_mex.hpp"
#else
#include "blasType.hpp"
#endif
#include "colMat.hpp"
#include "arrVec.hpp"
#include "skewMat.hpp"

class HouseholderMatrix
{
    typedef HouseholderMatrix SELF_TYPE;
    typedef ColMat<REAL_TYPE> MATX_TYPE;
    typedef ArrVec<REAL_TYPE> VECT_TYPE;

public:
    INTE_TYPE d;
    MATX_TYPE MatV;
    MATX_TYPE MatW;
    MATX_TYPE MatH;
    VECT_TYPE VecT;

    HouseholderMatrix() : d(0), MatV(ColMat<REAL_TYPE>()), MatW(ColMat<REAL_TYPE>()), MatH(ColMat<REAL_TYPE>()), VecT(ArrVec<REAL_TYPE>()) {};
    HouseholderMatrix(INTE_TYPE dim) : d(dim), MatV(ColMat<REAL_TYPE>(dim, dim - 1)), MatW(ColMat<REAL_TYPE>(dim, dim - 1)), MatH(ColMat<REAL_TYPE>(dim, dim)), VecT(ArrVec<REAL_TYPE>(dim - 1)) {};
    void copy(const SELF_TYPE &src)
    {
        // assert(this->d == src.d);
        MatV.copy(src.MatV);
        MatW.copy(src.MatW);
        MatH.copy(src.MatH);
        VecT.copy(src.VecT);
    }
    HouseholderMatrix(const SELF_TYPE &src) : HouseholderMatrix(src.d) { this->copy(src); };
    void swap(SELF_TYPE &rhs)
    {
        std::swap(this->d, rhs.d);
        this->MatV.swap(rhs.MatV);
        this->MatW.swap(rhs.MatW);
        this->MatH.swap(rhs.MatH);
        this->VecT.swap(rhs.VecT);
    };
    HouseholderMatrix &operator=(const SELF_TYPE &src)
    {
        SELF_TYPE temp(src);
        swap(temp);
        return (*this);
    };
    ~HouseholderMatrix() {};

    void Assign(REAL_TYPE *MatM, INTE_TYPE ldm)
    {
        MatH.Assign(MatM, ldm);
        MatV.Zero();
        MatW.Zero();
    };
    void Assign(ColMat<REAL_TYPE> MatM) { Assign(MatM.v, MatM.r); };

    void Copyto(SELF_TYPE &des) { des.copy(*this); };

    void Skew_Assign(REAL_TYPE *MatM, INTE_TYPE ldm)
    {
        MatH.Assign(MatM, ldm);
        for (auto ind = 0; ind < d; ind++)
            cblas_daxpby(d, -0.5, MatM + ind * ldm, 1, 0.5, MatH.v + ind, d);
        MatV.Zero();
        MatW.Zero();
    };
    void Skew_Assign(ColMat<REAL_TYPE> MatM) { Skew_Assign(MatM.v, MatM.r); };

    void Symm_Assign(REAL_TYPE *MatM, INTE_TYPE ldm)
    {
        MatH.Assign(MatM, ldm);
        for (auto ind = 0; ind < d; ind++)
            cblas_daxpby(d, 0.5, MatM + ind * ldm, 1, 0.5, MatH.v + ind, d);
        MatV.Zero();
        MatW.Zero();
    };
    void Symm_Assign(ColMat<REAL_TYPE> MatM) { Symm_Assign(MatM.v, MatM.r); };

    INTE_TYPE _vw_blk_size(INTE_TYPE n)
    {
        if (n <= 5)
            return 1;
        else if (n <= 12)
            return n - 4;
        else if (n <= 100)
            return 10;
        else
            return 60;
    }

    REAL_TYPE col_elimin(REAL_TYPE *_vec, REAL_TYPE *_col, INTE_TYPE n);
    // REAL_TYPE col_elimin(REAL_TYPE *_vec, INTE_TYPE n);

    void _SkewSymm_Hessenberg_VT(REAL_TYPE *MatV, INTE_TYPE ldv, REAL_TYPE *MatS, INTE_TYPE lds, INTE_TYPE col_num);
    void _SkewSymm_Hessenberg_VT(INTE_TYPE col_num) { _SkewSymm_Hessenberg_VT(MatV.v, d, MatH.v, d, col_num); };

    void SkewSymm_TriDiag_BLAS2()
    {
        _SkewSymm_Hessenberg_VT(MatV.v, d, MatH.v, d, d - 1);
        MatV.v[d * (d - 2) - 1] = 1.0;
        // VecT.v[d - 2] = 0.0;
    };
    void SkewSymm_TriDiag_BLAS2(REAL_TYPE *MatS, INTE_TYPE lds) { _SkewSymm_Hessenberg_VT(MatV.v, d, MatS, lds, d - 1); };

    void _SkewSymm_Hessenberg_VW(REAL_TYPE *MatV, INTE_TYPE ldv, REAL_TYPE *MatW, INTE_TYPE ldw, REAL_TYPE *MatS, INTE_TYPE lds, REAL_TYPE *VecTau, INTE_TYPE reduced_dim, INTE_TYPE col_num);
    void _SkewSymm_Hessenberg_VW(INTE_TYPE col_beg, INTE_TYPE col_num) { _SkewSymm_Hessenberg_VW(MatV.v, d, MatW.v, d, MatH.v, d, VecT.v, d - col_beg - 1, col_num); };
    void _SkewSymm_VW_Update(REAL_TYPE *MatS, INTE_TYPE lds, INTE_TYPE col_beg, INTE_TYPE col_num)
    {
        INTE_TYPE left_over_size = d - col_beg - col_num;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, left_over_size, left_over_size, col_num,
                    1.0, MatV.v + col_beg * (d + 1) + col_num, d, MatW.v + col_beg * (d + 1) + col_num, d,
                    1.0, MatS, lds);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, left_over_size, left_over_size, col_num,
                    -1.0, MatW.v + col_beg * (d + 1) + col_num, d, MatV.v + col_beg * (d + 1) + col_num, d,
                    1.0, MatS, lds);
    }

    void SkewSymm_TriDiag_BLAS3(REAL_TYPE *MatV, INTE_TYPE ldv, REAL_TYPE *MatW, INTE_TYPE ldw, REAL_TYPE *MatS, INTE_TYPE lds, REAL_TYPE *VecTau, INTE_TYPE nb)
    {
        // std::printf("SkewSymm_TriDiag_BLAS3 Call:\t Size:\t%lli, Block:\t%lli.\n", d, nb);

        INTE_TYPE col_beg = 0;

        REAL_TYPE *ptrV = MatV, *ptrW = MatW, *ptrS = MatS, *ptrT = VecTau;

        // View_ColMat<REAL_TYPE>(ptrS, lds, d, d).printf("Active MatS that is being Diagonalized:\n");

        while (col_beg + nb <= d - 2)
        {
            // View_ColMat<REAL_TYPE>(ptrS, lds, d - col_beg, d - col_beg).printf("Actived sub-matrix:\n");

            _SkewSymm_Hessenberg_VW(ptrV, ldv, ptrW, ldw, ptrS, lds, ptrT, d - col_beg - 1, nb);

            // View_ColMat<REAL_TYPE>(ptrS, lds, d - col_beg, d - col_beg).printf("Actived sub-matrix after elimination:\n");

            ptrV += nb * (ldv + 1);
            ptrW += nb * (ldw + 1);
            ptrS += nb * (lds + 1);
            ptrT += nb;
            _SkewSymm_VW_Update(ptrS, lds, col_beg, nb);

            col_beg += nb;
        }
        nb = d - 1 - col_beg;
        if (nb > 0)
        {
            _SkewSymm_Hessenberg_VW(ptrV, ldv, ptrW, ldw, ptrS, lds, ptrT, d - col_beg - 1, nb);
            ptrS += nb * (lds + 1);
            _SkewSymm_VW_Update(ptrS, lds, col_beg, nb);
        }
    }
    void SkewSymm_TriDiag_BLAS3(REAL_TYPE *MatS, INTE_TYPE lds, INTE_TYPE nb) { SkewSymm_TriDiag_BLAS3(MatV.v, d, MatW.v, d, MatS, lds, VecT.v, nb); };
    void SkewSymm_TriDiag_BLAS3(REAL_TYPE *MatS, INTE_TYPE lds) { SkewSymm_TriDiag_BLAS3(MatV.v, d, MatW.v, d, MatS, lds, VecT.v, _vw_blk_size(d)); };
    void SkewSymm_TriDiag_BLAS3(INTE_TYPE nb) { SkewSymm_TriDiag_BLAS3(MatV.v, d, MatW.v, d, MatH.v, d, VecT.v, nb); };
    void SkewSymm_TriDiag_BLAS3() { SkewSymm_TriDiag_BLAS3(MatV.v, d, MatW.v, d, MatH.v, d, VecT.v, _vw_blk_size(d)); };

    void Symm_TriDiag(REAL_TYPE *diag, REAL_TYPE *offd) { LAPACKE_dsytrd(LAPACK_COL_MAJOR, 'L', d, MatH.v, d, diag, offd, VecT.v); };

    void Symm_TriDiag() { LAPACKE_dsytrd(LAPACK_COL_MAJOR, 'L', d, MatH.v, d, MatV.v, MatV.v + d, VecT.v); };

    void Compute_DenseH(REAL_TYPE *H, INTE_TYPE ldh)
    {
        MatH.Zero();
        memcpy(MatH.v, MatV.v, sizeof(REAL_TYPE) * MatV.r * MatV.c);
        if (H != MatH.v)
            MatH.Copyto(H, ldh);
        LAPACKE_dorgtr(LAPACK_COL_MAJOR, 'L', d, H, ldh, VecT.v);
    };
    void Compute_DenseH() { Compute_DenseH(MatH.v, d); };

    void LeftAction(REAL_TYPE *MatX, INTE_TYPE colx, INTE_TYPE ldx, CHAR_TYPE trans) { LAPACKE_dormtr(LAPACK_COL_MAJOR, 'L', 'L', trans, d, colx, MatV.v, d, VecT.v, MatX, ldx); };
    void RightAction(REAL_TYPE *MatX, INTE_TYPE rowx, INTE_TYPE ldx, CHAR_TYPE trans) { LAPACKE_dormtr(LAPACK_COL_MAJOR, 'R', 'L', trans, rowx, d, MatV.v, d, VecT.v, MatX, ldx); };

    void Action(REAL_TYPE *MatX, INTE_TYPE colx, INTE_TYPE ldx) { LeftAction(MatX, colx, ldx, 'N'); };
    void Action(ColMat<REAL_TYPE> &MatX) { Action(MatX.v, MatX.c, MatX.r); };
    void Congruence(REAL_TYPE *MatX, INTE_TYPE ldx, CHAR_TYPE trans)
    {
        // Compute QXQ^T when trans = 'N' and Q^TXQ when trans = 'T';
        if (trans == 'N')
        {
            this->LeftAction(MatX, d, ldx, 'N');
            this->RightAction(MatX, d, ldx, 'T');
        }
        else if (trans == 'T')
        {
            this->LeftAction(MatX, d, ldx, 'T');
            this->RightAction(MatX, d, ldx, 'N');
        }
        else
        {
            // std::printf("Error! TRANS flag not supported!\n");
            throw(1);
        }
    };
    void Congruence(REAL_TYPE *MatX, INTE_TYPE ldx)
    {
        // Compute QXQ^T;
        this->LeftAction(MatX, d, ldx, 'N');
        this->RightAction(MatX, d, ldx, 'T');
    }

    void Copy_OffDiagonal(REAL_TYPE *offdiag)
    {
        auto ptrH = MatH.v + 1;
        for (INTE_TYPE ind = 0; ind < d - 1; ind++, ptrH += (d + 1), offdiag++)
            offdiag[0] = ptrH[0];
    }

    // void printf()
    // {
    //     std::printf("Householder Reflectors in compact V form:\n");
    //     MatV.printf();
    //     std::printf("Vector of tau values in Householder Reflectors:\n");
    //     for (INTE_TYPE ind = 0; ind < d - 1; ind++)
    //         std::printf("%f\n", VecT.v[ind]);
    //     // std::printf("Householder reflectors in compact W form (may be not computed):\n");
    //     // MatW.printf();
    //     std::printf("Householder reflectors in Matrix I - WV' form (may be not computed):\n");
    //     MatH.printf();
    // };
};

// class HouseholderReflector : public ColMat<REAL_TYPE>
// {
//     // This object stores the (n x n)-HouseholderReflector R = I - tau * v * v' in forms of scalar s and real vector v.
//     // By construction, v has leading 1 which reserve space for . The tau and v are aligned as (n)-vector denoted as (s, v_2, ..., v_n).
//     // It accepts multiple reflectors, stored in each columns.
// public:
//     REAL_TYPE *tau;
//     INTE_TYPE d;
//     INTE_TYPE *ind_os; // Offset of s, v in the colmun major matrix.
//     INTE_TYPE *el_len; // Length of v in the column major matrix.
//     BOOL_TYPE bl;      // Indicate if the v is stored in the lower triangular order, typically from a QR decomposition.

//     HouseholderReflector() : ColMat<REAL_TYPE>(), tau(nullptr), d(0), ind_os(nullptr), el_len(nullptr), bl(false) {};
//     HouseholderReflector(INTE_TYPE r, INTE_TYPE c) : ColMat<REAL_TYPE>(r, c), tau(new REAL_TYPE[c]), d(r), ind_os(nullptr), el_len(nullptr), bl(false) {};
//     void copy(const HouseholderReflector &src)
//     {
//         // ColMat<REAL_TYPE>::copy(src);
//         memcpy(this->v, src.v, sizeof(REAL_TYPE) * src.r * src.c);
//         memcpy(this->tau, src.tau, sizeof(REAL_TYPE) * c);
//         d = src.d;
//         this->bl = src.bl;
//         if (src.ind_os)
//             memcpy(this->ind_os, src.ind_os, sizeof(INTE_TYPE) * d);
//         if (src.el_len)
//             memcpy(this->el_len, src.el_len, sizeof(INTE_TYPE) * d);
//     };
//     HouseholderReflector(const HouseholderReflector &src) : HouseholderReflector(src.r, src.c)
//     {
//         if (src.ind_os)
//             this->ind_os = new INTE_TYPE[d];
//         if (src.el_len)
//             this->el_len = new INTE_TYPE[d];
//         this->copy(src);
//     }
//     void swap(HouseholderReflector &src)
//     {
//         ColMat<REAL_TYPE>::swap(src);
//         using std::swap;
//         swap(this->tau, src.tau);
//         swap(this->d, src.d);
//         swap(this->ind_os, src.ind_os);
//         swap(this->el_len, src.el_len);
//         swap(this->bl, src.bl);
//     }
//     HouseholderReflector &operator=(const HouseholderReflector &rhs)
//     {
//         HouseholderReflector temp(rhs);
//         swap(temp);
//         return (*this);
//     };

//     ~HouseholderReflector()
//     {
//         delete[] tau;
//         if (ind_os)
//             delete[] ind_os;
//         if (el_len)
//             delete[] el_len;
//     }

//     REAL_TYPE col_elimin(REAL_TYPE *vec, INTE_TYPE col, INTE_TYPE offset, INTE_TYPE length);
//     void Action(CHAR_TYPE side, CHAR_TYPE trans, REAL_TYPE *MatQ, INTE_TYPE ldq);

//     void printf()
//     {
//         std::printf("Elementary reflectors (vector):\n");
//         ColMat<REAL_TYPE>::printf();
//         std::printf("Elememtary reflectors (tau):\n");
//         for (INTE_TYPE ind = 0; ind < c; ind++)
//             std::printf("%f\t", tau[ind]);
//         std::printf("\n");
//     }
// };