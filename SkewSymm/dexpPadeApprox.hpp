#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

#include "colMajMat.hpp"
#include "skewSymmMat.hpp"
#include "pivotedLUFac.hpp"
#include "expmpa.hpp"
#include "dgemm.hpp"

// This implementation does not consider the scaling and squaring process, at this point.

#define _DEXPM_PADE_APPROXIMANT_PARANUM 13

class dexpPadeApproximant
{
public:
    matExpPadeApproximant *empa;
    ColMat<REAL_TYPE> **Ms;
    // Pade Approximant of A = S^{1/2^s} stored in epms->Ms, labelled as:   A, A2, A4, A6, W1, W2, W, Z1, Z2, U, V
    // The epa object does not belong to dexpPadeSeriesPara and should never be released in dexpPadeApproximant.
    ColMat<REAL_TYPE> **Rs;
    // Necessary intermediate matrices R, R^2, ..., R^{2^s} = Q = exp(S)
    // Rs is initialized on demand, with the maximum capacity specified in _EXPM_PADE_MAXIMUM_SCALING_ORDER.
    // This object belong the dexpPadeApproximant and should be managed by it.
    // In case of reusing the allocated Rs spaces, the dimension of the matrices are assumed to be exactly the same.
    ColMat<REAL_TYPE> **Ws;
    // Perturb-related Matrices, labelled as:                               M, M2, M4, M6, LW1, LW2, LZ1, LZ2, LW, LU, LV, R, L
    // Ws belongs to dexpPadeApproximant and should be managed by it. However, there is no need to deep copy Ws
    // as they are workspace assumed to be carrying trash unless inside a specific routine.

    INTE_TYPE d;

    dexpPadeApproximant()
    {
        empa = nullptr;
        Ms = nullptr;
        Rs = nullptr;
        Ws = nullptr;
        d = 0;
    };
    dexpPadeApproximant(INTE_TYPE dim)
    {
        d = dim;
        empa = nullptr;
        Ms = nullptr;
        Ws = new ColMat<REAL_TYPE> *[_DEXPM_PADE_APPROXIMANT_PARANUM];
        for (INTE_TYPE i = 0; i < _DEXPM_PADE_APPROXIMANT_PARANUM; i++)
            Ws[i] = new ColMat<REAL_TYPE>(d, d);

        Rs = new ColMat<REAL_TYPE> *[_EXPM_PADE_MAXIMUM_SCALING_ORDER];
        for (INTE_TYPE i = 0; i < _EXPM_PADE_MAXIMUM_SCALING_ORDER; i++)
            Rs[i] = nullptr;
    };
    dexpPadeApproximant(const dexpPadeApproximant &src) : dexpPadeApproximant(src.d)
    {
        this->d = src.d;
        this->empa = src.empa;
        this->Ms = src.Ms;

        // The workspace Ws are not copied.

        for (INTE_TYPE i = 0; i < _EXPM_PADE_MAXIMUM_SCALING_ORDER; i++)
            if (src.Rs[i])
            {
                this->Rs[i] = new ColMat<REAL_TYPE>(d, d);
                memcpy(this->Rs[i]->v, src.Rs[i]->v, sizeof(REAL_TYPE) * d * d);
            }
    };
    void swap(dexpPadeApproximant &src)
    {
        using std::swap;
        swap(this->empa, src.empa);
        swap(this->Ms, src.Ms);
        swap(this->Rs, src.Rs);
        swap(this->Ws, src.Ws);
    };
    dexpPadeApproximant &operator=(dexpPadeApproximant &rhs)
    {
        dexpPadeApproximant temp(rhs);
        swap(temp);
        return (*this);
    };
    ~dexpPadeApproximant()
    {
        // if (Ms)
        // {
        //     for (INTE_TYPE i = 0; i < _DEXPM_PADE_APPROXIMANT_PARANUM; i++)
        //         delete Ms[i];
        //     delete[] Ms;
        // }
        empa = nullptr;
        Ms = nullptr;
        if (Rs)
        {
            for (INTE_TYPE i = 0; i < _EXPM_PADE_MAXIMUM_SCALING_ORDER; i++)
                if (Rs[i])
                    delete Rs[i];
            delete[] Rs;
        }
        Rs = nullptr;

        if (Ws)
        {
            for (INTE_TYPE i = 0; i < _DEXPM_PADE_APPROXIMANT_PARANUM; i++)
                delete Ws[i];
            delete[] Ws;
        }
        Ws = nullptr;
    };

    void setup_empa(matExpPadeApproximant *_empa)
    {
        this->empa = _empa;
        this->Ms = empa->Ms;
    };
    void setup_rs();
    void setup_perturb(REAL_TYPE *M, INTE_TYPE ldm)
    {
        // assert(empa != nullptr)
        REAL_TYPE *src_ptr = M;
        REAL_TYPE *des_ptr = Ws[0]->v;
        for (INTE_TYPE col_ind = 0; col_ind < d; col_ind++, src_ptr += ldm, des_ptr += d)
            memcpy(des_ptr, src_ptr, sizeof(REAL_TYPE) * d);

        if (empa->s != 0)
            scal((1.0 / pow(2, empa->s)), *Ws[0]);
    };
    void setup_perturb(const ColMat<REAL_TYPE> &M) { setup_perturb(M.v, M.r); };

    void _dexp_core();
    void _dexp_recover(ColMat<REAL_TYPE> &N);
    void _dexp_recover(ColMat<REAL_TYPE> &N, ColMat<REAL_TYPE> &Q);

    void dexp(ColMat<REAL_TYPE> &N, const ColMat<REAL_TYPE> &M, matExpPadeApproximant *_empa)
    {
        setup_empa(_empa);
        setup_rs();
        setup_perturb(M);
        _dexp_core();
        _dexp_recover(N);
    };

    void dexp(ColMat<REAL_TYPE> &N, const ColMat<REAL_TYPE> &M)
    {
        setup_perturb(M);
        _dexp_core();
        _dexp_recover(N);
    };

    void dexp(ColMat<REAL_TYPE> &N)
    {
        _dexp_core();
        _dexp_recover(N);
    };

    void dexp(ColMat<REAL_TYPE> &N, ColMat<REAL_TYPE> &Q)
    {
        _dexp_core();
        _dexp_recover(N, Q);
    }

    friend void dexpPadeApprox(ColMat<REAL_TYPE> &N, ColMat<REAL_TYPE> &M, dexpPadeApproximant &Para) { Para.dexp(N, M); };
};