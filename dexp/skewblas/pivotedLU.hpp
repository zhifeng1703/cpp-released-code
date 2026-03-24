#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

#include "colMat.hpp"
#include "arrVec.hpp"

class PivotedLUFactor
{
    typedef PivotedLUFactor SELF_TYPE;
    typedef ColMat<REAL_TYPE> MATX_TYPE;
    typedef ArrVec<INTE_TYPE> VECT_TYPE;

public:
    MATX_TYPE M;
    VECT_TYPE P;
    INTE_TYPE r;
    INTE_TYPE c;
    INTE_TYPE s;
    BOOL_TYPE f;

    PivotedLUFactor() : M(ColMat<REAL_TYPE>()), P(ArrVec<INTE_TYPE>()) {};
    PivotedLUFactor(INTE_TYPE row, INTE_TYPE col) : M(ColMat<REAL_TYPE>(row, col)), P(ArrVec<INTE_TYPE>(row)), r(row), c(col), s(0), f(false) {};
    PivotedLUFactor(INTE_TYPE dim) : PivotedLUFactor(dim, dim) {};
    void copy(const SELF_TYPE &src)
    {
        this->M.copy(src.M);
        this->P.copy(src.P);
        this->r = src.r;
        this->c = src.c;
        this->s = src.s;
        this->f = src.f;
    }

    PivotedLUFactor(const View_ColMat<REAL_TYPE> &ViewM) : PivotedLUFactor(ViewM.r, ViewM.c) { Assign(ViewM.v, ViewM.ld); };
    PivotedLUFactor(const MATX_TYPE &MatM) : PivotedLUFactor(MatM.r, MatM.c) { Assign(MatM.v, MatM.r); };
    PivotedLUFactor(const SELF_TYPE &src) : PivotedLUFactor(src.r, src.c) { this->copy(src); };
    void swap(SELF_TYPE &src)
    {
        this->M.swap(src.M);
        this->P.swap(src.P);
        using std::swap;
        swap(this->r, src.r);
        swap(this->c, src.c);
        swap(this->s, src.s);
        swap(this->f, src.f);
    };
    PivotedLUFactor &operator=(const SELF_TYPE &rhs)
    {
        SELF_TYPE temp(rhs);
        swap(temp);
        return (*this);
    };
    ~PivotedLUFactor() {};

    void Assign(REAL_TYPE *MatA, INTE_TYPE lda)
    {
        M.Assign(MatA, lda);
        f = false;
    };
    void Assign(const View_ColMat<REAL_TYPE> &ViewA) { Assign(ViewA.v, ViewA.ld); };
    void Assign(const MATX_TYPE &MatA) { Assign(MatA.v, MatA.r); };

    INTE_TYPE Factor()
    {
        if (!f)
            s = LAPACKE_dgetrf(CblasColMajor, r, c, M.v, r, P.v);
        f = true;
        return s;
    }

    void Solve(REAL_TYPE *X, INTE_TYPE ldx, INTE_TYPE colx) { LAPACKE_dgetrs(CblasColMajor, 'N', r, colx, M.v, r, P.v, X, ldx); };
    void Solve(const View_ColMat<REAL_TYPE> ViewX) { Solve(ViewX.v, ViewX.ld, ViewX.c); };
};