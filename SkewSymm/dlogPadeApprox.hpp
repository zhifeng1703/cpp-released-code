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
#include "daxpy.hpp"
#include "dgemm.hpp"
// #include "dlaqtr.hpp"
#include "sylvester.hpp"

#define _DLOGM_PADE_APPROXIMANT_ORDER 13
#define _LOGM_PADE_MAXIMUM_SCALING_ORDER 8

class dlogPadeApproximant
{
    typedef SchurAngularFactor FACTOR;

public:
    FACTOR *saf;
    REAL_TYPE *a;
    REAL_TYPE *v;

    ColMat<REAL_TYPE> **Ts;
    ColMat<REAL_TYPE> **Es;
    ColMat<REAL_TYPE> **Ys;
    ColMat<REAL_TYPE> **Ws;

    REAL_TYPE *sys;

    INTE_TYPE d;
    INTE_TYPE s;
    INTE_TYPE asize;

    dlogPadeApproximant()
    {
        saf = nullptr;
        a = nullptr;
        v = nullptr;
        Ts = nullptr;
        Es = nullptr;
        Ys = nullptr;
        Ws = nullptr;
        sys = nullptr;
        d = 0;
        s = 0;
        asize = 0;
    };
    dlogPadeApproximant(INTE_TYPE dim)
    {
        d = dim;
        s = 0;
        asize = 0;

        Ts = new ColMat<REAL_TYPE> *[_LOGM_PADE_MAXIMUM_SCALING_ORDER];
        for (INTE_TYPE i = 0; i < _LOGM_PADE_MAXIMUM_SCALING_ORDER; i++)
            Ts[i] = new ColMat<REAL_TYPE>(d, d);
        Es = new ColMat<REAL_TYPE> *[_LOGM_PADE_MAXIMUM_SCALING_ORDER + 1];
        for (INTE_TYPE i = 0; i < _LOGM_PADE_MAXIMUM_SCALING_ORDER + 1; i++)
            Es[i] = new ColMat<REAL_TYPE>(d, d);
        Ys = new ColMat<REAL_TYPE> *[_DLOGM_PADE_APPROXIMANT_ORDER];
        for (INTE_TYPE i = 0; i < _DLOGM_PADE_APPROXIMANT_ORDER; i++)
            Ys[i] = new ColMat<REAL_TYPE>(d, d);
        Ws = new ColMat<REAL_TYPE> *[4];
        for (INTE_TYPE i = 0; i < 4; i++)
            Ws[i] = new ColMat<REAL_TYPE>(d, d);

        sys = new REAL_TYPE[24];
    };
    dlogPadeApproximant(const dlogPadeApproximant &src) : dlogPadeApproximant(src.d){};
    void swap(dlogPadeApproximant &src)
    {
        using std::swap;
        swap(this->saf, src.saf);
        swap(this->a, src.a);
        swap(this->v, src.v);

        swap(this->Ts, src.Ts);
        swap(this->Es, src.Es);
        swap(this->Ys, src.Ys);
        swap(this->Ws, src.Ws);

        swap(this->sys, src.sys);

        swap(this->d, src.d);
        swap(this->s, src.s);
        swap(this->asize, src.asize);
    };
    dlogPadeApproximant &operator=(dlogPadeApproximant &rhs)
    {
        dlogPadeApproximant temp(rhs);
        swap(temp);
        return (*this);
    };
    ~dlogPadeApproximant()
    {
        if (Ts)
        {
            for (INTE_TYPE i = 0; i < _LOGM_PADE_MAXIMUM_SCALING_ORDER; i++)
                if (Ts[i])
                    delete Ts[i];
            delete[] Ts;
        }
        Ts = nullptr;

        if (Es)
        {
            for (INTE_TYPE i = 0; i < _LOGM_PADE_MAXIMUM_SCALING_ORDER + 1; i++)
                if (Es[i])
                    delete Es[i];
            delete[] Es;
        }
        Es = nullptr;

        if (Ys)
        {
            for (INTE_TYPE i = 0; i < _DLOGM_PADE_APPROXIMANT_ORDER; i++)
                if (Ys[i])
                    delete Ys[i];
            delete[] Ys;
        }
        Ys = nullptr;

        if (Ws)
        {
            for (INTE_TYPE i = 0; i < 4; i++)
                delete Ws[i];
            delete[] Ws;
        }
        Ws = nullptr;

        if (sys)
            delete[] sys;
        sys = nullptr;
    };

    void setupSAF(FACTOR *factor)
    {
        saf = factor;
        v = saf->v;
        a = saf->a;
        asize = saf->asize;
    }

    void setup_scale()
    {
        s = 0;
        REAL_TYPE scale = 1.0;
        REAL_TYPE angle = 1.0;
        while (s < 3)
        {
            memset(Ts[s]->v, 0, sizeof(REAL_TYPE) * d * d);
            for (INTE_TYPE i = 0; i < asize; i++)
            {
                angle = a[i] / scale;
                Ts[s]->var(2 * i, 2 * i) = cos(angle);
                Ts[s]->var(2 * i + 1, 2 * i) = sin(angle);
                Ts[s]->var(2 * i, 2 * i + 1) = -sin(angle);
                Ts[s]->var(2 * i + 1, 2 * i + 1) = cos(angle);
            }
            s += 1;
            scale *= 2.0;
        }
    }

    void setup_perturb(REAL_TYPE *M, INTE_TYPE ldm)
    {
        // assert(empa != nullptr)

        Ws[0]->assign(M, ldm);
        my_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, d, d, d, 1.0, v, d, Ws[0]->v, d, 0.0, Ws[1]->v, d);
        my_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, 1.0, Ws[1]->v, d, v, d, 0.0, Ws[2]->v, d);

        Es[0]->assign(*Ws[2]);
    };

    void setup_perturb(REAL_TYPE *M, INTE_TYPE ldm, REAL_TYPE *Q, INTE_TYPE ldQ)
    {
        // assert(empa != nullptr)

        Ws[0]->assign(M, ldm);
        my_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, 1.0, Q, ldQ, Ws[0]->v, d, 0.0, Ws[1]->v, d);
        my_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, d, d, d, 1.0, v, d, Ws[1]->v, d, 0.0, Ws[2]->v, d);
        my_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, 1.0, Ws[2]->v, d, v, d, 0.0, Ws[3]->v, d);

        Es[0]->assign(*Ws[3]);
    };

    void setup_perturb(const ColMat<REAL_TYPE> &M) { setup_perturb(M.v, M.r); };

    void setup_perturb(const ColMat<REAL_TYPE> &M, const ColMat<REAL_TYPE> &Q) { setup_perturb(M.v, M.r, Q.v, Q.r); };

    void solve_Es()
    {
        // solve_Es_SYL();
        solve_Es_BLK();
    }

    void solve_Es_SYL()
    {
        ColMat<REAL_TYPE> **currE = Es;
        ColMat<REAL_TYPE> **currT = Ts;
        ColMat<REAL_TYPE> **nextE = Es + 1;
        for (INTE_TYPE i = 0; i < s; i++, currE++, currT++, nextE++)
        {
            (*nextE)->assign(**currE);
            symsyl(1, **nextE, **currT);
        }
    }

    void solve_Es_BLK()
    {
        ColMat<REAL_TYPE> **currE = Es;
        ColMat<REAL_TYPE> **currT = Ts;
        ColMat<REAL_TYPE> **nextE = Es + 1;

        REAL_TYPE sys_a, sys_b, sys_c, sys_d;
        REAL_TYPE *curr_blk;
        for (INTE_TYPE i = 0; i < s; i++, currE++, currT++, nextE++)
        {
            (*nextE)->assign(**currE);
            for (INTE_TYPE col_blk_ind = 0; col_blk_ind < asize; col_blk_ind++)
                for (INTE_TYPE row_blk_ind = 0; row_blk_ind < asize; row_blk_ind++)
                {
                    sys_a = (*currT)->ele(2 * row_blk_ind, 2 * row_blk_ind);
                    sys_b = (*currT)->ele(2 * row_blk_ind + 1, 2 * row_blk_ind);
                    sys_c = (*currT)->ele(2 * col_blk_ind, 2 * col_blk_ind);
                    sys_d = (*currT)->ele(2 * col_blk_ind + 1, 2 * col_blk_ind);
                    *(sys + 0) = sys_a + sys_c;
                    *(sys + 1) = -sys_b;
                    *(sys + 2) = sys_d;
                    *(sys + 3) = 0;
                    *(sys + 4) = sys_b;
                    *(sys + 5) = sys_a + sys_c;
                    *(sys + 6) = 0;
                    *(sys + 7) = sys_d;
                    *(sys + 8) = -sys_d;
                    *(sys + 9) = 0;
                    *(sys + 10) = sys_a + sys_c;
                    *(sys + 11) = -sys_b;
                    *(sys + 12) = 0;
                    *(sys + 13) = -sys_d;
                    *(sys + 14) = sys_b;
                    *(sys + 15) = sys_a + sys_c;

                    curr_blk = (*nextE)->pvar(2 * row_blk_ind, 2 * col_blk_ind);
                    *(sys + 16) = *(curr_blk + 0);
                    *(sys + 17) = *(curr_blk + 1);
                    *(sys + 18) = *(curr_blk + d);
                    *(sys + 19) = *(curr_blk + d + 1);

                    // LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', 4, 4, 1, sys, 4, sys + 16, 4);
                    cblas_dgemv(CblasColMajor, CblasNoTrans, 4, 4, 1.0, sys, 4, sys + 16, 1, 0.0, sys + 20, 1);

                    *(curr_blk + 0) = *(sys + 20);
                    *(curr_blk + 1) = *(sys + 21);
                    *(curr_blk + d) = *(sys + 22);
                    *(curr_blk + d + 1) = *(sys + 23);
                }
        }
    }

    void _dlog_core()
    {
        Ws[0]->assign(*Ts[s - 1]);
        for (INTE_TYPE ind = 0; ind < d; ind++)
            Ws[0]->var(ind, ind) -= 1.0;
        // get R in Ws[0];
        Ws[1]->assign(*Ws[0]);
        scal(0.1, *Ws[1]);
        for (INTE_TYPE ind = 0; ind < d; ind++)
            Ws[1]->var(ind, ind) += 1.0;
        // get scaled R in Ws[1];

        ColMat<REAL_TYPE> *currY = nullptr;
        for (INTE_TYPE i = 0; i < _LOGM_PADE_MAXIMUM_SCALING_ORDER; i++)
        {
            currY = Ys[i];
            currY->assign(*Es[s]);
            // my_dlaqtr_real(false, d, Ws[1]->v, d, 1.0, currY->v, Ws[2]->v);
            // my_dlaqtr_real(true, d, Ws[1]->v, d, 1.0, currY->v, Ws[2]->v);
            // LAPACKE_dtrtrs(CblasColMajor, 'U', 'N', 'N', d, d, Ws[1]->v, d, currY->v, d);
            // LAPACKE_dtrtrs(CblasRowMajor, 'U', 'T', 'N', d, d, Ws[1]->v, d, currY->v, d);
        }
    }
    void _dlog_recover(ColMat<REAL_TYPE> &N)
    {

        memset(Ws[0]->v, 0, sizeof(REAL_TYPE) * d * d);
        ColMat<REAL_TYPE> **currY = Ys;
        for (INTE_TYPE i = 0; i < _LOGM_PADE_MAXIMUM_SCALING_ORDER; i++, currY++)
            my_daxpy(d * d, 1.0, (*currY)->v, Ws[0]->v);
        my_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, d, d, d, 1.0, v, d, Ws[0]->v, d, 0.0, Ws[1]->v, d);
        my_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, d, d, d, 1.0, Ws[1]->v, d, v, d, 0.0, N.v, d);
        my_dscal(d * d, pow(2.0, s), N.v);
    };

    void dlog(ColMat<REAL_TYPE> &N, const ColMat<REAL_TYPE> &M)
    {
        setup_perturb(M);
        setup_scale();
        solve_Es();
        _dlog_core();
        _dlog_recover(N);
    };

    void dlog(ColMat<REAL_TYPE> &N, const ColMat<REAL_TYPE> &M, const ColMat<REAL_TYPE> &Q)
    {
        setup_perturb(M, Q);
        setup_scale();
        solve_Es();
        _dlog_core();
        _dlog_recover(N);
    };

    // void dlog(ColMat<REAL_TYPE> &N, const ColMat<REAL_TYPE> &M)
    //{
    //     setup_perturb(M);
    //     solve_Es();
    //     _dlog_core();
    //     _dlog_recover(N);
    // };
};