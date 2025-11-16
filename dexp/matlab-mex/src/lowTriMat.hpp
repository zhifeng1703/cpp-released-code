#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>

#include "colMat.hpp"
#include "lowTra.hpp"

template <class ELEM_TYPE>
class LowTriMat : public ColMat<ELEM_TYPE>
{
    typedef LowTriMat<ELEM_TYPE> SELF_TYPE;

public:
    LowerTraversal tra; // Traversal of the lower triangular part of the matrix. By default it use the column traversal that takes no extra space
    ELEM_TYPE *lv;      // Extra vector for the LowerTri part. By default it is not initialized.
    INTE_TYPE d;        // Dimension of the matrix
    INTE_TYPE msize;    // # of the entries in the matrix vector
    INTE_TYPE lsize;    // # of the entries in the lowerTri vector
    BOOL_TYPE sl;       // Flag that indicate if the diagonal is stored.

    LowTriMat() : ColMat<ELEM_TYPE>(), tra(LowerTraversal(0, COL_DIAG)), lv(nullptr), d(0), msize(0), lsize(0), sl(false) {};
    LowTriMat(INTE_TYPE dim) : ColMat<ELEM_TYPE>(dim, dim), tra(LowerTraversal(dim, COL_DIAG)), lv(nullptr), d(dim), msize(dim * dim), lsize(0), sl(false) {};
    LowTriMat(INTE_TYPE dim, LOWTRA_TYPE lt, BOOL_TYPE strict_lower) : ColMat<ELEM_TYPE>(dim, dim), tra(LowerTraversal(dim, lt)), lv(nullptr), d(dim), msize(dim * dim), lsize(0), sl(strict_lower) {};
    void init_low_vec()
    {
        lsize = sl ? ((d - 1) * d / 2) : ((d + 1) * d / 2);
        if (!lv)
            lv = new ELEM_TYPE[((d + 1) * d / 2)];
    };
    void copy(const SELF_TYPE &src)
    {
        ColMat<ELEM_TYPE>::copy(src);
        this->tra.copy(src.tra);
        this->d = src.d;
        this->msize = src.msize;
        this->lsize = src.lsize;
        this->sl = src.sl;
        if (src.lv)
        {
            this->init_low_vec();
            memcpy(this->lv, src.lv, sizeof(ELEM_TYPE) * lsize);
        }
    }
    LowTriMat(const SELF_TYPE &src) : LowTriMat(src.d)
    {
        this->copy(src);
    };
    void swap(SELF_TYPE &rhs)
    {
        ColMat<ELEM_TYPE>::swap(rhs);
        using std::swap;
        swap(this->tra, rhs.tra);
        swap(this->lv, rhs.lv);
        swap(this->d, rhs.d);
        swap(this->msize, rhs.msize);
        swap(this->lsize, rhs.lsize);
        swap(this->sl, rhs.sl);
    };
    LowTriMat<ELEM_TYPE> &operator=(const SELF_TYPE &rhs)
    {
        SELF_TYPE temp(rhs);
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

    void mat2vec(const LowerTraversal &lt) { lt.mat2vec(lv, this->v, d); };
    void mat2vec() { mat2vec(tra); };

    void vec2low(const LowerTraversal &lt) { lt.vec2low(this->v, d, lv); };
    void vec2low() { vec2low(tra); };

    void vec2upp(const LowerTraversal &lt) { lt.vec2upp(this->v, d, lv); };
    void vec2upp() { vec2upp(tra); };

    // void vec2low(INTE_TYPE beg, INTE_TYPE end)
    // {
    //     for (auto i = beg; i <= end; i++)
    //         ColMat<ELEM_TYPE>::v[tra->v2l[i]] = lv[i];
    // };

    void Zero()
    {
        if (sl)
        {
            ELEM_TYPE *col_ptr = this->v + 1;
            for (auto col_ind = 0; col_ind < d - 1; col_ind++, col_ptr += d + 1)
                memset(col_ptr, 0, sizeof(ELEM_TYPE) * (d - 1 - col_ind));
        }
        else
        {
            ELEM_TYPE *col_ptr = this->v;
            for (auto col_ind = 0; col_ind < d; col_ind++, col_ptr += d + 1)
                memset(col_ptr, 0, sizeof(ELEM_TYPE) * (d - col_ind));
        }
    }

    template <class random_engine>
    void Rand(random_engine &engine, ELEM_TYPE a, ELEM_TYPE b)
    {
        std::uniform_real_distribution<ELEM_TYPE> distribution(a, b);
        if (sl)
        {
            ELEM_TYPE *col_ptr = this->v + 1;
            for (auto col_ind = 0; col_ind < d - 1; col_ind++, col_ptr += d + 1)
                for (auto row_ind = 0; row_ind < d - col_ind - 1; row_ind++)
                    col_ptr[row_ind] = distribution(engine);
        }
        else
        {
            ELEM_TYPE *col_ptr = this->v;
            for (auto col_ind = 0; col_ind < d; col_ind++, col_ptr += d + 1)
                for (auto row_ind = 0; row_ind < d - col_ind; row_ind++)
                    col_ptr[row_ind] = distribution(engine);
        }
    }
};
