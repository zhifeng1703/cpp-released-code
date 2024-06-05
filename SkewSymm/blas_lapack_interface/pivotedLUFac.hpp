#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

#include "type_convention.hpp"
#include "dgetrf.hpp"
#include "dgetrs.hpp"

class PivotedLUFactor
{
public:
    ColMat<REAL_TYPE> M;
    INTE_TYPE *P;
    INTE_TYPE r;
    INTE_TYPE c;
    INTE_TYPE s;
    BOOL_TYPE f;

    PivotedLUFactor()
    {
        M = ColMat<REAL_TYPE>();
        P = nullptr;
        r = 0;
        c = 0;
        s = 0;
        f = false;
    };
    PivotedLUFactor(INTE_TYPE row, INTE_TYPE col) : M(ColMat<REAL_TYPE>(row, col)), P(new INTE_TYPE[row]), r(row), c(col), s(0), f(false){};
    PivotedLUFactor(INTE_TYPE dim) : PivotedLUFactor(dim, dim){};
    PivotedLUFactor(const View_ColMat<REAL_TYPE> &ViewM) : PivotedLUFactor(ViewM.r, ViewM.c)
    {
        for (INTE_TYPE col_ind = 0; col_ind < c; col_ind++)
            memcpy(M[col_ind], ViewM.v + (col_ind * ViewM.ld), sizeof(REAL_TYPE) * r);
    };
    PivotedLUFactor(const PivotedLUFactor &src) : PivotedLUFactor(src.r, src.c)
    {
        memcpy(this->M.v, src.M.v, sizeof(REAL_TYPE) * r * c);
        if (src.P)
            memcpy(this->P, src.P, sizeof(INTE_TYPE) * r);
    };
    void swap(PivotedLUFactor &src)
    {
        this->M.swap(src.M);
        using std::swap;
        swap(this->P, src.P);
        swap(this->r, src.r);
        swap(this->c, src.c);
        swap(this->s, src.s);
        swap(this->f, src.f);
    };
    PivotedLUFactor &operator=(const PivotedLUFactor &rhs)
    {
        PivotedLUFactor temp(rhs);
        swap(temp);
        return (*this);
    };
    ~PivotedLUFactor()
    {
        if (P)
            delete[] P;
        P = nullptr;
    };

    void copy(const PivotedLUFactor &src)
    {
        memcpy(this->M.v, src.M.v, sizeof(REAL_TYPE) * r * c);
        if (src.P)
            memcpy(this->P, src.P, sizeof(INTE_TYPE) * r);
    }

    void release_factor()
    {
        if (P)
            delete[] P;
        P = nullptr;
    };

    void setupMat(REAL_TYPE *A, INTE_TYPE lda)
    {
        for (INTE_TYPE col_ind = 0; col_ind < c; col_ind++)
            memcpy(M[col_ind], A + col_ind * lda, sizeof(REAL_TYPE) * r);
        f = false;
    };
    void setupMat(REAL_TYPE *A)
    {
        memcpy(M.v, A, sizeof(REAL_TYPE) * r * c);
        f = false;
    };
    void setupMat(const View_ColMat<REAL_TYPE> &ViewA) { setupMat(ViewA.v, ViewA.ld); };
    void setupMat(const ColMat<REAL_TYPE> &MatA) { setupMat(MatA.v); };

    INTE_TYPE computePLU()
    {
        if (!f)
            s = my_dgetrf(r, c, M.v, r, P);
        return s;
    }

    friend void solverPLU(REAL_TYPE *X, INTE_TYPE ldX, INTE_TYPE colX, const PivotedLUFactor &plu);
    friend void solverPLU(const View_ColMat<REAL_TYPE> ViewX, const PivotedLUFactor &plu);
};