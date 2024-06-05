#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>

#include "type_convention.hpp"
#include "lowTriMat.hpp"
#include "schAngFac.hpp"
#include "spectralFac.hpp"
#include "dgees.hpp"

class SkewSymmMat : public LowTriMat<REAL_TYPE>
{
    typedef SchurAngularFactor FACTOR;

public:
    FACTOR saf; // The factor is not a necessary object, it only allocate when the decomposition routine is called.

    SkewSymmMat() : LowTriMat<REAL_TYPE>(), saf(FACTOR()){};
    SkewSymmMat(INTE_TYPE dim) : LowTriMat<REAL_TYPE>(dim), saf(FACTOR()){};
    void copy(const SkewSymmMat &src)
    {
        LowTriMat<REAL_TYPE>::copy(src);
        if (!src.saf.is_empty())
            this->initial_saf();
        saf.copy(src.saf);
    };
    SkewSymmMat(const SkewSymmMat &src) : SkewSymmMat(src.d) { copy(src); };
    void swap(SkewSymmMat &src)
    {
        this->LowTriMat::swap(src);
        this->saf.swap(src.saf);
    };
    SkewSymmMat &operator=(const SkewSymmMat &src)
    {
        SkewSymmMat temp(src);
        // this->SkewSymmMat::swap(temp);
        swap(temp);
        return (*this);
    };
    ~SkewSymmMat(){};

    void setdiag()
    {
        for (INTE_TYPE i = 0; i < d; i++)
            ColMat<REAL_TYPE>::v[i + i * d] = 0.0;
    };

    void vec2upp()
    {
        for (INTE_TYPE i = 0; i < lsize; i++)
            ColMat<REAL_TYPE>::v[tra->v2u[i]] = -lv[i];
    };

    void vec2mat()
    {
        vec2low();
        vec2upp();
        setdiag();
    };

    void initial_saf(REAL_TYPE *work)
    {
        if (saf.is_empty())
            saf = FACTOR(d, work);
    };
    void initial_saf()
    {
        if (saf.is_empty())
            saf = FACTOR(d);
    };

    void SchurAngular_SkewSymm();

    friend REAL_TYPE *Generate_Workspace(const SkewSymmMat *S)
    {
        return new REAL_TYPE[S->d * (S->d * 2)];
    }

    friend REAL_TYPE *Generate_Workspace(const SkewSymmMat *S, INTE_TYPE dim)
    {
        return new REAL_TYPE[dim * (dim + 2)];
    }

    friend void SchurAngular_SkewSymm(SchurAngularFactor &saf, const SkewSymmMat &MatS);

    friend void SchurAngular_SkewSymm_Spectral(SchurAngularFactor &saf, const SkewSymmMat &MatS);

    friend void Spectral_SkewSymm(SpectralFactor &spf, const SkewSymmMat &MatS);
};
