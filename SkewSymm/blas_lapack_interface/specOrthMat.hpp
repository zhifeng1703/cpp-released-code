#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

#include "type_convention.hpp"
#include "schAngFac.hpp"
#include "dgees.hpp"

class SpecOrthMat : public ColMat<REAL_TYPE>
{
    typedef SchurAngularFactor FACTOR;

public:
    FACTOR saf;
    INTE_TYPE d;

    SpecOrthMat() : ColMat<REAL_TYPE>(), saf(FACTOR()), d(0){};
    SpecOrthMat(INTE_TYPE dim) : ColMat<REAL_TYPE>(dim, dim), saf(FACTOR()), d(dim){};
    void copy(const SpecOrthMat &src)
    {
        ColMat<REAL_TYPE>::copy(src);
        if (!src.saf.is_empty())
        {
            this->initial_saf();
            saf.copy(src.saf);
        }
    }
    SpecOrthMat(const SpecOrthMat &src) : ColMat<REAL_TYPE>(src.d, src.d) { copy(src); };
    void swap(SpecOrthMat &src)
    {
        ColMat<REAL_TYPE>::swap(src);
        this->saf.swap(src.saf);
        std::swap(this->d, src.d);
    };
    SpecOrthMat &operator=(const SpecOrthMat &rhs)
    {
        SpecOrthMat temp(rhs);
        this->swap(temp);
        return (*this);
    };

    void initial_saf()
    {
        if (saf.is_empty())
            saf = FACTOR(d);
    };
    void initial_saf(REAL_TYPE *work)
    {
        if (saf.is_empty())
            saf = FACTOR(d, work);
    };

    void SchurAngular_SpecOrth();
    friend void SchurAngular_SpecOrth(FACTOR &saf, const SpecOrthMat &MatS);
};