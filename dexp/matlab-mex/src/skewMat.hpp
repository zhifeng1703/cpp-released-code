#pragma once

#include <cassert>
#include <cstring>

#ifdef MATLAB_MEX_BUILD
#include "blasType_mex.hpp"
#else
#include "blasType.hpp"
#endif
#include "lowTriMat.hpp"

class SkewSymmMat : public LowTriMat<REAL_TYPE>
{
public:
    SkewSymmMat() : LowTriMat<REAL_TYPE>() {};
    SkewSymmMat(INTE_TYPE dim) : LowTriMat<REAL_TYPE>(dim, COL_OFFD, true) {};
    void copy(const SkewSymmMat &src) { LowTriMat<REAL_TYPE>::copy(src); };
    SkewSymmMat(const SkewSymmMat &src) : SkewSymmMat(src.d) { copy(src); };
    void swap(SkewSymmMat &src) { this->LowTriMat::swap(src); };
    SkewSymmMat &operator=(const SkewSymmMat &src)
    {
        SkewSymmMat temp(src);
        swap(temp);
        return (*this);
    };
    ~SkewSymmMat() {};

    void setdiag()
    {
        for (auto i = 0; i < d; i++)
            ColMat<REAL_TYPE>::v[i + i * d] = 0.0;
    };

    void vec2upp(const LowerTraversal &lt)
    {
        cblas_dscal(lsize, -1, lv, 1);
        LowTriMat::vec2upp(lt);
        cblas_dscal(lsize, -1, lv, 1);
    };

    void vec2mat(const LowerTraversal &lt)
    {
        vec2low(lt);
        vec2upp(lt);
        setdiag();
    };

    void low2upp(const LowerTraversal &lt)
    {
        mat2vec(lt);
        cblas_dscal(lsize, -1, lv, 1);
        LowTriMat::vec2upp(lt);
        cblas_dscal(lsize, -1, lv, 1);
    };

    void vec2upp()
    {
        cblas_dscal(lsize, -1, lv, 1);
        LowTriMat::vec2upp(tra);
        cblas_dscal(lsize, -1, lv, 1);
    };

    void low2upp()
    {
        if (!lv)
            this->init_low_vec();
        mat2vec(tra);
        cblas_dscal(lsize, -1, lv, 1);
        LowTriMat::vec2upp(tra);
        cblas_dscal(lsize, -1, lv, 1);
    };

    void vec2mat() { vec2mat(tra); };
};
