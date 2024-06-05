#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cmath>

#include "type_convention.hpp"
#include "dgees.hpp"

class SchurCanonicalFactor
{
public:
    ColMat<REAL_TYPE> svec;
    ColMat<REAL_TYPE> uppt;
    REAL_TYPE *valr;
    REAL_TYPE *vali;
    INTE_TYPE d;
    INTE_TYPE nzsize;

    SchurCanonicalFactor() : svec(ColMat<REAL_TYPE>()), uppt(ColMat<REAL_TYPE>()), valr(nullptr), vali(nullptr), d(0), nzsize(0){};
    SchurCanonicalFactor(INTE_TYPE dim) : svec(ColMat<REAL_TYPE>(dim, dim)), uppt(ColMat<REAL_TYPE>(dim, dim)), valr(new REAL_TYPE[dim]), vali(new REAL_TYPE[dim]), d(dim), nzsize(0){
        svec.fast_col_access();
    };
    SchurCanonicalFactor(const SchurCanonicalFactor &src) : SchurCanonicalFactor(src.d)
    {
        this->svec.copy(src.svec);
        this->uppt.copy(src.uppt);
        memcpy(this->valr, src.valr, sizeof(REAL_TYPE) * d);
        memcpy(this->vali, src.vali, sizeof(REAL_TYPE) * d);
        this->nzsize = src.nzsize;
    };
    void swap(SchurCanonicalFactor &rhs)
    {
        this->svec.swap(rhs.svec);
        this->uppt.swap(rhs.uppt);
        using std::swap;
        swap(this->valr, rhs.valr);
        swap(this->vali, rhs.vali);
        swap(this->d, rhs.d);
        swap(this->nzsize, rhs.nzsize);
    };
    SchurCanonicalFactor &operator=(const SchurCanonicalFactor &rhs)
    {
        SchurCanonicalFactor temp(rhs);
        swap(temp);
        return (*this);
    };
    ~SchurCanonicalFactor()
    {
        if (valr)
            delete[] valr;
        valr = nullptr;
        if (vali)
            delete[] vali;
        vali = nullptr;
    };

    void setup_decompose_matrix(REAL_TYPE* M, INTE_TYPE ldm)
    {
        REAL_TYPE *src = M;
        REAL_TYPE *des = uppt.v;
        for (INTE_TYPE col_ind = 0; col_ind < d; col_ind++, des += d, src += ldm)
            memcpy(des, src, sizeof(REAL_TYPE) * d);
    };
    void setup_decompose_matrix(REAL_TYPE* M)
    {
        memcpy(uppt.v, M, sizeof(REAL_TYPE) * d * d);
    };
    void setup_decompose_matrix(const View_ColMat<REAL_TYPE> &ViewM) { setup_decompose_matrix(ViewM.v, ViewM.ld); };
    void setup_decompose_matrix(const ColMat<REAL_TYPE> &MatM) { setup_decompose_matrix(MatM.v); };

    void factor(){
        // Assuming uppt is assigned with the matrix that needs to be factored.
        my_dgees('V', 'S', _dgees_select_nonzero, d, uppt.v, d, &nzsize, valr, vali, svec.v, d);
    };
    void factor(REAL_TYPE* M, INTE_TYPE ldm) {
        setup_decompose_matrix(M, ldm);
        factor();
    };
    void factor(REAL_TYPE* M) {
        setup_decompose_matrix(M);
        factor();
    };    
    void factor(const View_ColMat<REAL_TYPE> &ViewM) {
        setup_decompose_matrix(ViewM);
        factor();
    };
    void factor(const ColMat<REAL_TYPE> &MatM) {
        setup_decompose_matrix(MatM);
        factor();
    };
};