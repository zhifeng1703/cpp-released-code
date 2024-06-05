#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>

#include "type_convention.hpp"
#include "colMajMat.hpp"
#include "lowerTraversal.hpp"

template <class type>
class LowTriMat : public ColMat<type>
{
public:
    LowerTraversal *tra; // Traversal of the lower triangular part of the matrix.
    type *lv;            // LowerTri vector
    INTE_TYPE d;         // Dimension of the matrix
    INTE_TYPE msize;     // # of the entries in the matrix vector
    INTE_TYPE lsize;     // # of the entries in the lowerTri vector
    BOOL_TYPE sl;        // Flag that indicate if the diagonal is stored.

    LowTriMat() : ColMat<type>(), tra(nullptr), lv(nullptr), d(0), msize(0), lsize(0), sl(false){};
    LowTriMat(INTE_TYPE dim) : ColMat<type>(dim, dim), tra(nullptr), lv(nullptr), d(dim), msize(dim * dim), lsize(0), sl(false){};
    void copy(const LowTriMat<type> &src)
    {
        ColMat<type>::copy(src);
        this->tra = src.tra;
        this->d = src.d;
        this->msize = src.msize;
        this->lsize = src.lsize;
        this->sl = src.sl;
        if (src.lv)
        {
            this->init_low_vec(this->sl);
            memcpy(this->lv, src.lv, sizeof(type) * (d * (d + 1) / 2));
        }
    }
    LowTriMat(const LowTriMat<type> &src) : LowTriMat(src.d)
    {
        this->copy(src);
    };
    void swap(LowTriMat<type> &rhs)
    {
        ColMat<type>::swap(rhs);
        using std::swap;
        swap(this->tra, rhs.tra);
        swap(this->lv, rhs.lv);
        swap(this->d, rhs.d);
        swap(this->msize, rhs.msize);
        swap(this->lsize, rhs.lsize);
        swap(this->sl, rhs.sl);
    };
    LowTriMat<type> &operator=(const LowTriMat<type> &rhs)
    {
        LowTriMat<type> temp(rhs);
        swap(temp);
        return (*this);
    };
    ~LowTriMat()
    {
        // v belongs to mat and will be released over there.
        if (lv)
            delete[] lv;
        lv = nullptr;
    };

    void mat2vec()
    {
        for (INTE_TYPE i = 0; i < lsize; i++)
            lv[i] = ColMat<type>::v[tra->v2l[i]];
    };

    void vec2low()
    {
        for (INTE_TYPE i = 0; i < lsize; i++)
            ColMat<type>::v[tra->v2l[i]] = lv[i];
    };

    void init_low_vec(BOOL_TYPE strict_lower)
    {
        sl = strict_lower;
        lsize = sl ? ((d - 1) * d) / 2 : ((d + 1) * d / 2);
        if (!lv)
            lv = new type[((d + 1) * d / 2)];
    };

    void set_low_vec(LowerTraversal *traversal)
    {
        tra = traversal;
        sl = tra->sl;
        lsize = sl ? ((d - 1) * d) / 2 : ((d + 1) * d / 2);
        if (!lv)
            lv = new type[((d + 1) * d / 2)];
        this->mat2vec();
    };
};
