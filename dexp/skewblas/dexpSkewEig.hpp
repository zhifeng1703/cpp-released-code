#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

#include "colMat.hpp"
#include "arrVec.hpp"
#include "skewSchFac.hpp"
#include "skewMat.hpp"

class dexpSkewEigenPara
{
    typedef dexpSkewEigenPara SELF_TYPE;
    typedef ArrVec<CMPX_TYPE> VECT_TYPE;
    typedef ColMat<CMPX_TYPE> MATX_TYPE;

public:
    // Native members owned by dexpSemiSimplePara:

    INTE_TYPE d;

    MATX_TYPE _forward_para;
    MATX_TYPE _inverse_para;
    MATX_TYPE EigVec;
    VECT_TYPE EigVal;

    MATX_TYPE MatX;
    MATX_TYPE MatY;
    MATX_TYPE Work;

    // CMPX_TYPE *forward;
    // CMPX_TYPE *inverse;
    // CMPX_TYPE *work;
    // INTE_TYPE d;

    dexpSkewEigenPara() {};
    dexpSkewEigenPara(INTE_TYPE dim) : _forward_para(MATX_TYPE(dim, dim)), _inverse_para(MATX_TYPE(dim, dim)),
                                       EigVec(MATX_TYPE(dim, dim)), EigVal(VECT_TYPE(dim)),
                                       MatX(MATX_TYPE(dim, dim)), MatY(MATX_TYPE(dim, dim)), Work(MATX_TYPE(dim, dim)), d(dim) {};
    void copy(const SELF_TYPE &src)
    {
        // assert(this->normal == src.normal)
        this->_forward_para.copy(src._forward_para);
        this->_inverse_para.copy(src._inverse_para);
        this->EigVec.copy(src.EigVec);
        this->EigVal.copy(src.EigVal);
    };
    void swap(SELF_TYPE &src)
    {
        this->_forward_para.swap(src._forward_para);
        this->_inverse_para.swap(src._inverse_para);
        this->EigVec.swap(src.EigVec);
        this->EigVal.swap(src.EigVal);
        this->MatX.swap(src.MatX);
        this->MatY.swap(src.MatY);
        this->Work.swap(src.Work);
        using std::swap;
        swap(this->d, src.d);
    };
    dexpSkewEigenPara(const SELF_TYPE &src) : dexpSkewEigenPara(src.d) { this->copy(src); };
    SELF_TYPE &operator=(const SELF_TYPE &src)
    {
        SELF_TYPE temp(src);
        swap(temp);
        return (*this);
    };
    ~dexpSkewEigenPara() {};

    void Assign(const ColMat<CMPX_TYPE> &vectors, const ArrVec<CMPX_TYPE> &values)
    {
        this->EigVec.Assign(vectors);
        this->EigVal.Assign(values.v);
    }

    void Assign(const SkewSchurFactor &ssf)
    {
        REAL_TYPE *Q = reinterpret_cast<REAL_TYPE *>(EigVec.v);
        REAL_TYPE *D = reinterpret_cast<REAL_TYPE *>(EigVal.v);

        memset(D, 0, sizeof(REAL_TYPE) * 2 * d);
        cblas_daxpy(ssf.a, -1.0, ssf.A.v, 1, D + 1, 4);
        cblas_dcopy(ssf.a, ssf.A.v, 1, D + 3, 4);

        memset(Q, 0, sizeof(REAL_TYPE) * 2 * d * d);
        INTE_TYPE m = d / 2;
        REAL_TYPE s = sqrt(2.0) / 2.0;
        for (auto pair_ind = 0; pair_ind < m; pair_ind++)
        {
            cblas_daxpy(d, s, ssf.R.v + (2 * pair_ind * d), 1, Q + (4 * pair_ind * d), 2);
            cblas_daxpy(d, s, ssf.R.v + (2 * pair_ind * d), 1, Q + ((2 * pair_ind + 1) * 2 * d), 2);
            cblas_daxpy(d, s, ssf.R.v + ((2 * pair_ind + 1) * d), 1, Q + (4 * pair_ind * d) + 1, 2);
            cblas_daxpy(d, -s, ssf.R.v + ((2 * pair_ind + 1) * d), 1, Q + ((2 * pair_ind + 1) * 2 * d) + 1, 2);
        }
        if (m + m != d)
            cblas_dcopy(d, ssf.R.v + (d - 1) * d, 1, Q + (d - 1) * 2 * d, 2);
    }

    void Parameter() { Parameter_Skew(); };
    void Parameter_Skew()
    {
        INTE_TYPE mat_ind = 0;
        CMPX_TYPE diff;
        auto _val = EigVal.v;
        auto _forward = _forward_para.v;
        auto _inverse = _inverse_para.v;
        for (INTE_TYPE col_ind = 0; col_ind < d; col_ind++)
            for (INTE_TYPE row_ind = 0; row_ind < d; row_ind++)
            {
                // diff =  val[col_ind] - val[row_ind];
                substrCMPX(diff, _val[col_ind], _val[row_ind]);
                if (normofCMPX(diff) < 1e-14)
                    assignCMPX(_forward[mat_ind], 1.0, 0.0);
                else
                {
                    // forward[mat_ind] = (exp(diff) - 1.0) / diff;
                    expm1divCMPX(_forward[mat_ind], diff);
                    // exponeCMPX(_forward[mat_ind], diff);
                    // _forward[mat_ind].real = _forward[mat_ind].real - 1.0;
                    // divideCMPX(_forward[mat_ind], _forward[mat_ind], diff);
                }
                // inverse[mat_ind] = 1.0 / forward[mat_ind];
                inversCMPX(_inverse[mat_ind], _forward[mat_ind]);
                mat_ind++;
            }
    };

    // void Parameter_Full()
    // {
    //     INTE_TYPE mat_ind = 0;
    //     CMPX_TYPE diff, num;
    //     REAL_TYPE diffnorm;
    //     auto _val = EigVal.v;
    //     auto _forward = _forward_para.v;
    //     auto _inverse = _inverse_para.v;
    //     auto eval_exp = ArrVec<CMPX_TYPE>(d);

    //     for (INTE_TYPE col_ind = 0; col_ind < d; col_ind++)
    //     {
    //         exponeCMPX(eval_exp[col_ind], _val[col_ind]);
    //         for (INTE_TYPE row_ind = 0; row_ind < d; row_ind++)
    //         {
    //             // diff =  val[col_ind] - val[row_ind];
    //             // REAL_TYPE scale = 1.0 + fabs(_val[col_ind].real) + fabs(_val[col_ind].imag) + fabs(_val[row_ind].real) + fabs(_val[row_ind].imag);
    //             // REAL_TYPE tol0 = 1e-14 * scale;
    //             // REAL_TYPE tol1 = 1e-8 * scale;
    //             substrCMPX(diff, _val[col_ind], _val[row_ind]);
    //             diffnorm = normofCMPX(diff);
    //             if (diffnorm < 1e-14)
    //                 assignCMPX(_forward[mat_ind], eval_exp[col_ind]);
    //             else if (diffnorm < 1e-8)
    //             {
    //                 CMPX_TYPE diff2;
    //                 multifCMPX(diff2, diff, diff); // diff^2

    //                 // term = diff/2 + diff^2/6
    //                 CMPX_TYPE t1, t2;
    //                 multifCMPX(t1, 0.5, diff);
    //                 multifCMPX(t2, 1.0 / 6.0, diff2);
    //                 additiCMPX(num, t1, t2);

    //                 assignCMPX(_forward[mat_ind], 1.0 + num.real, num.imag);             // res = 1 + term
    //                 multifCMPX(_forward[mat_ind], eval_exp[row_ind], _forward[mat_ind]); // e^y * res
    //             }
    //             else
    //             {
    //                 // forward[mat_ind] = (exp(diff) - 1.0) / diff;
    //                 // exponeCMPX(_forward[mat_ind], _val[col_ind]);
    //                 // exponeCMPX(_inverse[mat_ind], _val[row_ind]);
    //                 // substrCMPX(_forward[mat_ind], _forward[mat_ind], _inverse[mat_ind]);
    //                 // divideCMPX(_forward[mat_ind], _forward[mat_ind], diff);
    //                 substrCMPX(num, eval_exp[col_ind], eval_exp[row_ind]);
    //                 divideCMPX(_forward[mat_ind], num, diff);
    //             }
    //             // inverse[mat_ind] = 1.0 / forward[mat_ind];
    //             inversCMPX(_inverse[mat_ind], _forward[mat_ind]);
    //             mat_ind++;
    //         }
    //     }
    // };

    void Parameter_Full()
    {
        INTE_TYPE mat_ind = 0;
        CMPX_TYPE diff, diff2, term;
        REAL_TYPE diffnorm;
        auto _val = EigVal.v;
        auto _forward = _forward_para.v;
        auto _inverse = _inverse_para.v;
        auto eval_exp = ArrVec<CMPX_TYPE>(d);

        for (INTE_TYPE col_ind = 0; col_ind < d; ++col_ind)
            exponeCMPX(eval_exp[col_ind], _val[col_ind]);

        for (INTE_TYPE col_ind = 0; col_ind < d; ++col_ind)
            for (INTE_TYPE row_ind = 0; row_ind < d; ++row_ind, ++mat_ind)
            {
                substrCMPX(diff, _val[col_ind], _val[row_ind]);
                diffnorm = normofCMPX(diff);

                if (diffnorm < 1e-14)
                {
                    // diagonal: e^{λ_i}
                    assignCMPX(_forward[mat_ind], eval_exp[col_ind]);
                }
                else if (diffnorm < 1e-8)
                {
                    // small δ expansion: e^{λ_j}(1 + δ/2 + δ²/6)
                    multifCMPX(diff2, diff, diff); // δ²
                    multifCMPX(term, 0.5, diff);   // δ/2
                    CMPX_TYPE term2;
                    multifCMPX(term2, 1.0 / 6.0, diff2);
                    additiCMPX(term, term, term2);      // δ/2 + δ²/6
                    additiCMPX(term, term, {1.0, 0.0}); // 1 + ...
                    multifCMPX(_forward[mat_ind], eval_exp[row_ind], term);
                }
                else
                {
                    // general case: (exp(λ_i) - exp(λ_j)) / (λ_i - λ_j)
                    CMPX_TYPE num;
                    substrCMPX(num, eval_exp[col_ind], eval_exp[row_ind]);
                    divideCMPX(_forward[mat_ind], num, diff);
                }

                // build inverse stably: (λ_i - λ_j)/(exp(λ_i)-exp(λ_j))
                inversCMPX(_inverse[mat_ind], _forward[mat_ind]);
            }
    }

    void _assign_complex(CMPX_TYPE *CmatX, INTE_TYPE ldcx, REAL_TYPE *X, INTE_TYPE ldx)
    {
        REAL_TYPE *realC = reinterpret_cast<REAL_TYPE *>(CmatX);
        REAL_TYPE *ptrX = X;

        for (INTE_TYPE col = 0; col < d; col++, ptrX += ldx, realC += 2 * ldcx)
            cblas_dcopy(d, ptrX, 1, realC, 2);
    };

    void _retrieve_real(REAL_TYPE *X, INTE_TYPE ldx, CMPX_TYPE *CmatX, INTE_TYPE ldcx)
    {
        REAL_TYPE *realC = reinterpret_cast<REAL_TYPE *>(CmatX);
        REAL_TYPE *ptrX = X;

        for (INTE_TYPE col = 0; col < d; col++, ptrX += ldx, realC += 2 * ldcx)
            cblas_dcopy(d, realC, 2, ptrX, 1);
    };

    void _retrieve_norm(REAL_TYPE *X, INTE_TYPE ldx, CMPX_TYPE *CmatX, INTE_TYPE ldcx)
    {
        REAL_TYPE *ptrX = X;
        CMPX_TYPE *ptrC = CmatX;

        for (INTE_TYPE col = 0; col < d; col++, ptrX += ldx, ptrC += ldcx)
            for (INTE_TYPE row = 0; row < d; row++)
                ptrX[row] = normofCMPX(ptrC[row]);
    };

    void Dexp(REAL_TYPE *RealY, INTE_TYPE ldy, REAL_TYPE *RealX, INTE_TYPE ldx)
    {
        CMPX_TYPE complex0, complex1;
        assignCMPX(complex0, 0.0, 0.0);
        assignCMPX(complex1, 1.0, 0.0);

        _assign_complex(MatX.v, d, RealX, ldx);

        cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, d, d, d, &complex1, EigVec.v, d, MatX.v, d, &complex0, Work.v, d); // work = R^{H} * X;
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &complex1, Work.v, d, EigVec.v, d, &complex0, MatX.v, d);   // X = work * R;

        for (auto ind = 0; ind < d * d; ind++)
            multifCMPX(MatY.v[ind], _forward_para.v[ind], MatX.v[ind]);

        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &complex1, EigVec.v, d, MatY.v, d, &complex0, Work.v, d);   // work = R * Y;
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, d, d, d, &complex1, Work.v, d, EigVec.v, d, &complex0, MatY.v, d); // Y = work * R^{H};

        _retrieve_real(RealY, ldy, MatY.v, d);
        //_retrieve_norm(RealY, ldy, MatY.v, d);
    }

    void DexpInv(REAL_TYPE *RealY, INTE_TYPE ldy, REAL_TYPE *RealX, INTE_TYPE ldx)
    {
        CMPX_TYPE complex0, complex1;
        assignCMPX(complex0, 0.0, 0.0);
        assignCMPX(complex1, 1.0, 0.0);

        _assign_complex(MatX.v, d, RealX, ldx);

        cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, d, d, d, &complex1, EigVec.v, d, MatX.v, d, &complex0, Work.v, d); // work = R^{H} * X;
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &complex1, Work.v, d, EigVec.v, d, &complex0, MatX.v, d);   // X = work * R;

        for (auto ind = 0; ind < d * d; ind++)
            multifCMPX(MatY.v[ind], _inverse_para.v[ind], MatX.v[ind]);

        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, &complex1, EigVec.v, d, MatY.v, d, &complex0, Work.v, d);   // work = R * Y;
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, d, d, d, &complex1, Work.v, d, EigVec.v, d, &complex0, MatY.v, d); // Y = work * R^{H};

        _retrieve_real(RealY, ldy, MatY.v, d);
        // _retrieve_norm(RealY, ldy, MatY.v, d);
    }
};
