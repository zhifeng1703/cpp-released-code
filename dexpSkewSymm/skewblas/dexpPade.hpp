#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

#include "colMat.hpp"
#include "skewMat.hpp"
#include "pivotedLU.hpp"
#include "expmPade.hpp"

class dexpPadeApprox : public expmPadeApprox
{
    typedef expmPadeApprox EXPM_TYPE;
    typedef dexpPadeApprox SELF_TYPE;
    typedef ColMat<REAL_TYPE> MATX_TYPE;
    typedef ArrVec<ColMat<REAL_TYPE> *> VECT_TYPE;

public:
    VECT_TYPE VecMatE;
    // Perturb-related Matrices, labelled as:                               M, M2, M4, M6, LW1, LW2, LZ1, LZ2, LW, LU, LV, R, L
    VECT_TYPE VecMatD;

    // There is no need to deep copy Ws
    // as they are workspace assumed to be carrying trash unless inside a specific routine.

    dexpPadeApprox() : EXPM_TYPE()
    {
        VecMatE = VECT_TYPE();
        VecMatD = VECT_TYPE();
    };
    dexpPadeApprox(INTE_TYPE dim) : EXPM_TYPE(dim)
    {
        VecMatE = VECT_TYPE(_EXPM_PADE_APPROX_MAX_PARANUM);
        VecMatD = VECT_TYPE(_EXPM_PADE_APPROX_MAX_SCALING + 1);
    };
    void copy(const SELF_TYPE &src)
    {
        EXPM_TYPE::copy(src);
        _check_and_initialize(this->VecMatE, src._vmsize, src.d);
        for (auto ind = 0; ind < _vmsize; ind++)
            VecMatE[ind]->copy(*src.VecMatE[ind]);

        _check_and_initialize(this->VecMatD, src._vqsize, src.d);
        for (auto ind = 0; ind < _vmsize; ind++)
            VecMatD[ind]->copy(*src.VecMatD[ind]);
    }
    dexpPadeApprox(const SELF_TYPE &src) : SELF_TYPE(src.d) { this->copy(src); };
    void swap(SELF_TYPE &src)
    {
        EXPM_TYPE::swap(src);
        this->VecMatE.swap(src.VecMatE);
        this->VecMatD.swap(src.VecMatD);
    };
    SELF_TYPE &operator=(SELF_TYPE &rhs)
    {
        SELF_TYPE temp(rhs);
        swap(temp);
        return (*this);
    };
    ~dexpPadeApprox()
    {
        _release_vec_mat(VecMatE);
        _release_vec_mat(VecMatD);
    };

    void _assign_perturb(REAL_TYPE *MatE, INTE_TYPE lde)
    {
        //  Assuming the Pade approximant is computed and available.
        _check_and_initialize(VecMatE, 1, d);
        VecMatE[0]->Assign(MatE, lde);
        if (s != 0)
            VecMatE[0]->Scale(1.0 / pow(2, s));
    }
    void _assign_perturb(const ColMat<REAL_TYPE> &MatE) { _assign_perturb(MatE.v, MatE.r); };
    void _assign_perturb(const View_ColMat<REAL_TYPE> &ViewE) { _assign_perturb(ViewE.v, ViewE.ld); };

    void _dexp_pade_low();
    void _dexp_pade_13();

    void _dexp_pade()
    {
        if (m < 13)
            _dexp_pade_low();
        else
            _dexp_pade_13();
    };

    void _dexp_square()
    {
        // Recursively do L_{i+1} = L_iR_i + R_iL_i

        _check_and_initialize(VecMatD, s + 1, d);

        VecMatD[0]->Assign(*VecMatE[_vmsize - 1]);
        for (auto ind = 0; ind < s; ind++)
        {
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, 1.0, VecMatD[ind]->v, d, VecMatQ[ind]->v, d, 0.0, VecMatD[ind + 1]->v, d);
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, 1.0, VecMatQ[ind]->v, d, VecMatD[ind]->v, d, 1.0, VecMatD[ind + 1]->v, d);
        }
    }

    void Dexp(REAL_TYPE *MatX, INTE_TYPE ldx)
    {
        // Assuming Expm have been computed and available.
        _assign_perturb(MatX, ldx);
        _dexp_pade();
        _dexp_square();
    };
    void Dexp(MATX_TYPE &MatA) { Dexp(MatA.v, MatA.r); };
    void Dexp(View_ColMat<REAL_TYPE> &ViewA) { Dexp(ViewA.v, ViewA.ld); };

    void GetDexp(REAL_TYPE *MatY, INTE_TYPE ldy) { VecMatD[s]->Copyto(MatY, ldy); };

    void Dexp(REAL_TYPE *MatY, INTE_TYPE ldy, REAL_TYPE *MatX, INTE_TYPE ldx)
    {
        Dexp(MatX, ldx);
        GetDexp(MatY, ldy);
    }
    void Dexp(MATX_TYPE &MatY, MATX_TYPE &MatX) { Dexp(MatY.v, MatY.r, MatX.v, MatX.r); };
    void Dexp(View_ColMat<REAL_TYPE> &ViewY, View_ColMat<REAL_TYPE> &ViewX) { Dexp(ViewY.v, ViewY.ld, ViewX.v, ViewX.ld); };
};