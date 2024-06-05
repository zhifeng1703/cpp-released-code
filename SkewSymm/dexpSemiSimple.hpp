#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

#include "schAngFac.hpp"
#include "spectralFac.hpp"
#include "skewSymmMat.hpp"
#include "zgemm.hpp"

class dexpSemiSimplePara
{
    typedef SpectralFactor FACTOR;

public:
    // Native members owned by dexpSemiSimplePara:
    CMPX_TYPE *forward;
    CMPX_TYPE *inverse;
    CMPX_TYPE *work;
    INTE_TYPE d;

    // Friend member not owne by dexpSemiSimplePara:
    FACTOR *spf;
    CMPX_TYPE *_spv;
    CMPX_TYPE *_inv;
    CMPX_TYPE *_val;
    BOOL_TYPE _normal;

    dexpSemiSimplePara()
    {
        forward = nullptr;
        inverse = nullptr;
        work = nullptr;
        spf = nullptr;
        _spv = nullptr;
        _inv = nullptr;
        _val = nullptr;
        _normal = false;
    };
    dexpSemiSimplePara(INTE_TYPE dim) : forward(new CMPX_TYPE[dim * dim]), inverse(new CMPX_TYPE[dim * dim]), work(new CMPX_TYPE[2 * dim * dim]), d(dim)
    {
        spf = nullptr;
        _spv = nullptr;
        _inv = nullptr;
        _val = nullptr;
        _normal = false;
    };
    dexpSemiSimplePara(FACTOR *specfac) : dexpSemiSimplePara(specfac->d)
    {
        spf = specfac;
        _spv = spf->spv.v;
        _inv = spf->inv.v;
        _val = spf->val;
        _normal = spf->normal;
    };
    dexpSemiSimplePara(const dexpSemiSimplePara &src) : dexpSemiSimplePara(src.spf)
    {
        memcpy(this->forward, src.forward, d * d * sizeof(CMPX_TYPE));
        memcpy(this->inverse, src.inverse, d * d * sizeof(CMPX_TYPE));
    }
    void swap(dexpSemiSimplePara &src)
    {
        using std::swap;
        swap(this->forward, src.forward);
        swap(this->inverse, src.inverse);
        swap(this->work, src.work);
        swap(this->d, src.d);
        swap(this->spf, src.spf);
        swap(this->_spv, src._spv);
        swap(this->_inv, src._inv);
        swap(this->_val, src._val);
        swap(this->_normal, src._normal);
    };
    dexpSemiSimplePara &operator=(const dexpSemiSimplePara &rhs)
    {
        dexpSemiSimplePara temp(rhs);
        swap(temp);
        return (*this);
    };
    ~dexpSemiSimplePara()
    {
        if (forward)
            delete[] forward;
        forward = nullptr;
        if (inverse)
            delete[] inverse;
        inverse = nullptr;
        if (work)
            delete[] work;
        work = nullptr;
    };

    void setSPF(FACTOR *factor)
    {
        spf = factor;
        _spv = spf->spv.v;
        _inv = spf->inv.v;
        _val = spf->val;
        _normal = spf->normal;
    };

    void setupPara()
    {
        INTE_TYPE mat_ind = 0;
        CMPX_TYPE diff;
        for (INTE_TYPE col_ind = 0; col_ind < d; col_ind++)
            for (INTE_TYPE row_ind = 0; row_ind < d; row_ind++)
            {
                // diff =  val[col_ind] - val[row_ind];
                substrCMPX(diff, _val[col_ind], _val[row_ind]);
                if (normofCMPX(diff) < 1e-14)
                    assignCMPX(forward[mat_ind], 1.0, 0.0);
                else
                {
                    // forward[mat_ind] = (exp(diff) - 1.0) / diff;

                    exponeCMPX(forward[mat_ind], diff);
                    forward[mat_ind].real = forward[mat_ind].real - 1.0;
                    divideCMPX(forward[mat_ind], forward[mat_ind], diff);
                }
                // inverse[mat_ind] = 1.0 / forward[mat_ind];
                inversCMPX(inverse[mat_ind], forward[mat_ind]);
                mat_ind++;
            }
    };

    void _assignComplex(CMPX_TYPE *CmatX, INTE_TYPE ldcx, REAL_TYPE *X, INTE_TYPE ldx)
    {
        // INTE_TYPE mat_ind = 0;
        for (INTE_TYPE col_ind = 0; col_ind < d; col_ind++)
            for (INTE_TYPE row_ind = 0; row_ind < d; row_ind++ /*, mat_ind++*/)
                assignCMPX(CmatX[row_ind + col_ind * ldcx], X[row_ind + col_ind * ldx], 0.0);
    };
    void _assignComplex(CMPX_TYPE *CmatX, INTE_TYPE ldcx, const SkewSymmMat &MatX) { _assignComplex(CmatX, ldcx, MatX.v, MatX.d); };

    void _retrieveReal(REAL_TYPE *X, INTE_TYPE ldx, CMPX_TYPE *CmatX, INTE_TYPE ldcx)
    {
        // INTE_TYPE mat_ind = 0;
        for (INTE_TYPE col_ind = 0; col_ind < d; col_ind++)
            for (INTE_TYPE row_ind = 0; row_ind < d; row_ind++ /*, mat_ind++*/)
                X[row_ind + col_ind * ldx] = CmatX[row_ind + col_ind * ldcx].real;
    };
    void _retrieveReal(SkewSymmMat &MatX, CMPX_TYPE *CmatX, INTE_TYPE ldcx)
    {
        _retrieveReal(MatX.v, MatX.d, CmatX, ldcx);
        MatX.mat2vec();
    };

    void _dexpSemiSimple_similarity(CMPX_TYPE *Y, INTE_TYPE ldy, CMPX_TYPE *X, INTE_TYPE ldx, CMPX_TYPE *workS, INTE_TYPE ldw, BOOL_TYPE forward_action) const
    {
        CMPX_TYPE complex0, complex1;
        assignCMPX(complex0, 0.0, 0.0);
        assignCMPX(complex1, 1.0, 0.0);
        if (forward_action) // Y = R^{-1} X R
        {
            if (_normal)
            {
                my_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, d, d, d, complex1, _spv, d, X, ldx, complex0, workS, ldw); // work = R^{H} * X;
                my_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, complex1, workS, ldw, _spv, d, complex0, Y, ldy);   // Y = work * R;
            }
            else
            {
                my_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, complex1, _inv, d, X, ldx, complex0, workS, ldw); // work = R^{-1} * X;
                my_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, complex1, workS, ldw, _spv, d, complex0, Y, ldy); // Y = work * R;
            }
        }
        else // Y = R X R^{-1}
        {
            if (_normal)
            {
                my_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, complex1, _spv, d, X, ldx, complex0, workS, ldw);   // work = R * X;
                my_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, d, d, d, complex1, workS, ldw, _spv, d, complex0, Y, ldy); // Y = work * R^{H};
            }
            else
            {
                my_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, complex1, _spv, d, X, ldx, complex0, workS, ldw); // work = R * X;
                my_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, complex1, workS, ldw, _inv, d, complex0, Y, ldy); // Y = work * R';
            }
        }
    }
    void _dexpSemiSimple_core(CMPX_TYPE *Y, INTE_TYPE ldy, CMPX_TYPE *X, INTE_TYPE ldx, BOOL_TYPE forward_action)
    {
        auto map = forward_action ? forward : inverse;
        INTE_TYPE map_ind = 0;
        for (INTE_TYPE col_ind = 0; col_ind < d; col_ind++)
            for (INTE_TYPE row_ind = 0; row_ind < d; row_ind++)
                multifCMPX(Y[row_ind + col_ind * ldy], map[map_ind++], X[row_ind + col_ind * ldx]);
        // Y[row_ind + col_ind * ldy] = map[map_ind++] * X[row_ind + col_ind * ldx]
    }

    friend void dexpSemiSimple(REAL_TYPE *Y, INTE_TYPE ldy, REAL_TYPE *X, INTE_TYPE ldx, dexpSemiSimplePara &Para, BOOL_TYPE forward_action);
    friend void dexpSemiSimple(const View_ColMat<REAL_TYPE> &MatY, const View_ColMat<REAL_TYPE> &MatX, dexpSemiSimplePara &Para, BOOL_TYPE forward_action);
    friend void dexpSemiSimple(SkewSymmMat &MatY, SkewSymmMat &MatX, dexpSemiSimplePara &Para, BOOL_TYPE forward_action);
};
