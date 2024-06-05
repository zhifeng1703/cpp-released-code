#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

#include "type_convention.hpp"
#include "schAngFac.hpp"

class SpectralFactor
{
public:
    ColMat<CMPX_TYPE> spv;
    ColMat<CMPX_TYPE> inv;
    CMPX_TYPE *val;
    REAL_TYPE *w;
    INTE_TYPE d;
    BOOL_TYPE normal;

    SpectralFactor()
    {
        spv = ColMat<CMPX_TYPE>();
        inv = ColMat<CMPX_TYPE>();
        val = nullptr;
        w = nullptr;
        d = 0;
        normal = false;
    };
    SpectralFactor(INTE_TYPE dim, BOOL_TYPE is_normal) : val(new CMPX_TYPE[dim]), d(dim), normal(is_normal), w(new REAL_TYPE[2 * dim - 1])
    {
        spv = ColMat<CMPX_TYPE>(d, d);
        if (normal)
            inv = ColMat<CMPX_TYPE>();
        else
            inv = ColMat<CMPX_TYPE>(d, d);
    };
    SpectralFactor(const SchurAngularFactor &saf);
    SpectralFactor(const SpectralFactor &src)
    {
        spv = ColMat<CMPX_TYPE>(src.spv);
        inv = ColMat<CMPX_TYPE>(src.inv);
        val = new CMPX_TYPE[src.d];
        d = src.d;
        normal = src.normal;
        memcpy(val, src.val, d * sizeof(CMPX_TYPE));
    };
    void swap(SpectralFactor &src)
    {
        spv.swap(src.spv);
        inv.swap(src.inv);
        using std::swap;
        swap(val, src.val);
        swap(d, src.d);
        swap(normal, src.normal);
    };
    SpectralFactor &operator=(const SpectralFactor &rhs)
    {
        SpectralFactor temp(rhs);
        swap(temp);
        return (*this);
    };
    ~SpectralFactor()
    {
        if (val)
            delete[] val;
        val = nullptr;
        if (w)
            delete[] w;
        w = nullptr;
    };

    void setSPF(const SchurAngularFactor &saf);
};