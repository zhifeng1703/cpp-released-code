#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <algorithm>
#include <random>

#ifdef MATLAB_MEX_BUILD
#include "blasType_mex.hpp"
#else
#include "blasType.hpp"
#endif

#define COLMAT_DEFAULT_PRINT_COUNT 10
#define COLMAT_DEFAULT_PRINT_FORMAT "%1.3f\t"
#define MATRIX_INFINITE_NORM -1

template <class ELEM_TYPE>
class View_ColMat
{
    // The View object NEVER owns the array it accesses. Any View object must be initialized with the array it accesses.
    // It is encouraged to create and release View objects on demand for one-time-use, i.e., it is not recommended to reuse it,
    // i.e., the copy constructor and the assignment constructor are not recommended, as they are called only in reuse scenarios.
    // Also note that ColMat that does not own the vector v, i.e., b = false, functions exactly as View_ColMat with ld = r.

    typedef View_ColMat<ELEM_TYPE> SELF_TYPE;

public:
    ELEM_TYPE *v;
    INTE_TYPE ld; // The leading order of the view, which is the row number of the underlying ColMat or the equivalent structure.
    INTE_TYPE r;
    INTE_TYPE c;

    View_ColMat() : v(nullptr), ld(0), r(0), c(0) {};
    View_ColMat(ELEM_TYPE *vec, INTE_TYPE leading_order, INTE_TYPE row, INTE_TYPE col) : v(vec), ld(leading_order), r(row), c(col) {};
    void copy(const SELF_TYPE &src)
    {
        // no deep copy should be executed for View ELEM_TYPE.
        this->v = src.v;
        this->r = src.r;
        this->c = src.c;
        this->ld = src.ld;
    };
    View_ColMat(const SELF_TYPE &src) : View_ColMat(src.v, src.ld, src.r, src.c) {};
    void swap(SELF_TYPE &rhs)
    {
        using std::swap;
        swap(this->v, rhs.v);
        swap(this->r, rhs.r);
        swap(this->c, rhs.c);
        swap(this->ld, rhs.ld);
    }
    View_ColMat<ELEM_TYPE> &operator=(const SELF_TYPE &rhs)
    {
        View_ColMat<ELEM_TYPE> temp(rhs);
        swap(temp);
        return (*this);
    }
    ~View_ColMat() {};

    void Assign(ELEM_TYPE *src, INTE_TYPE lds, INTE_TYPE rs, INTE_TYPE cs)
    {

        if ((lds == ld) && (rs == lds) && (r == ld))
            memcpy(v, src, sizeof(ELEM_TYPE) * rs * cs);
        else
        {
            ELEM_TYPE *des_ptr = v;
            ELEM_TYPE *src_ptr = src;
            for (auto col_ind = 0; col_ind < cs; col_ind++, des_ptr += this->ld, src_ptr += lds)
                memcpy(des_ptr, src_ptr, sizeof(ELEM_TYPE) * rs);
        }
    }
    void Assign(const SELF_TYPE &src) { Assign(src.v, src.ld, src.r, src.c); };
    void Assign(ELEM_TYPE *src, INTE_TYPE lds) { Assign(src, lds, r, c); };
    void Assign(ELEM_TYPE *src) { Assign(src, r, r, c); };

    void Copyto(ELEM_TYPE *des, INTE_TYPE ldd, INTE_TYPE rd, INTE_TYPE cd)
    {
        if ((ldd == ld) && (rd == ldd) && (r == ld))
            memcpy(des, this->v, sizeof(ELEM_TYPE) * rd * cd);
        else
        {
            ELEM_TYPE *des_col = des;
            ELEM_TYPE *src_col = this->v;
            for (auto col_ind = 0; col_ind < cd; col_ind++, src_col += this->ld, des_col += ldd)
                memcpy(des_col, src_col, sizeof(ELEM_TYPE) * rd);
        }
    };
    void Copyto(const SELF_TYPE &des) { Copyto(des.v, des.ld, des.r, des.c); };
    void Copyto(ELEM_TYPE *des, INTE_TYPE ldd) { Copyto(des, ldd, r, c); };
    void Copyto(ELEM_TYPE *des) { Copyto(des, r, r, c); };

    ELEM_TYPE &operator()(INTE_TYPE i, INTE_TYPE j) { return *(v + i + j * ld); }; // return writable i,j entry
    ELEM_TYPE *operator[](INTE_TYPE i) { return v + i * ld; };                     // return the column vector

    ELEM_TYPE *ptr(INTE_TYPE i) { return v + i * ld; };
    ELEM_TYPE *ptr() { return v; };

    ELEM_TYPE &cstele(INTE_TYPE i, INTE_TYPE j) const { return *(v + i + j * ld); };
    ELEM_TYPE *cstptr(INTE_TYPE i, INTE_TYPE j) const { return v + i + j * ld; };

    void Reset(ELEM_TYPE *vec, INTE_TYPE leading_order, INTE_TYPE row, INTE_TYPE col)
    {
        // This routine basically functions as the re-assignment to lhs = View_ColMat<ELEM_TYPE>(vec, leading_order, row, col)
        // but it avoids the multiple constructor, thanks to the simple data structure.
        v = vec;
        ld = leading_order;
        r = row;
        c = col;
    };

    // It is very important to note that the (deep) copy operator of the View object is NOT called in the copy constructor
    // as this object MUST be an shallow object that does not own or operate any data unless specifically instructed.

    // void fprintf(FILE *of, const char *prefix, INTE_TYPE row_limit, INTE_TYPE col_limit, const char *format)
    // {
    //     using std::fprintf;

    //     INTE_TYPE row = (row_limit < r) ? row_limit : r;
    //     INTE_TYPE col = (col_limit < c) ? col_limit : c;

    //     fprintf(of, "%s", prefix);
    //     fprintf(of, "[\n");
    //     for (auto row_ind = 0; row_ind < row; row_ind++, (col == c ? fprintf(of, "%s", ";\n") : fprintf(of, "%s", "...\n")))
    //         for (auto col_ind = 0; col_ind < col; col_ind++)
    //             fprintf(of, format, ptr(col_ind)[row_ind]);
    //     if (row == r)
    //     {
    //         for (auto col_ind = 0; col_ind < col; col_ind++)
    //             fprintf(of, "%s", "...\t");
    //         fprintf(of, "%s", "\n");
    //     }
    //     fprintf(of, "%s", "]\n\n");
    // };

    // void fprintf(FILE *of, const char *prefix, INTE_TYPE row_limit, INTE_TYPE col_limit) { fprintf(of, prefix, row_limit, col_limit, COLMAT_DEFAULT_PRINT_FORMAT); };
    // void fprintf(FILE *of, const char *prefix, const char *format) { fprintf(of, prefix, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_COUNT, format); };
    // void fprintf(FILE *of, const char *prefix) { fprintf(of, prefix, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_FORMAT); };
    // void fprintf(FILE *of, INTE_TYPE row_limit, INTE_TYPE col_limit) { fprintf(of, "", row_limit, col_limit, COLMAT_DEFAULT_PRINT_FORMAT); };
    // void fprintf(FILE *of) { fprintf(of, "", COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_FORMAT); };

    // void printf(const char *prefix, INTE_TYPE row_limit, INTE_TYPE col_limit, const char *format)
    // {
    //     using std::printf;

    //     INTE_TYPE row = (row_limit < r) ? row_limit : r;
    //     INTE_TYPE col = (col_limit < c) ? col_limit : c;

    //     printf("%s", prefix);
    //     printf("[\n");
    //     for (auto row_ind = 0; row_ind < row; row_ind++, (col == c ? printf(";\n") : printf("...\n")))
    //         for (auto col_ind = 0; col_ind < col; col_ind++)
    //             printf(format, ptr(col_ind)[row_ind]);
    //     if (row != r)
    //     {
    //         for (auto col_ind = 0; col_ind < col; col_ind++)
    //             printf("...\t");
    //         printf("\n");
    //     }
    //     printf("]\n\n");
    // };

    // void printf(const char *prefix, INTE_TYPE row_limit, INTE_TYPE col_limit) { printf(prefix, row_limit, col_limit, COLMAT_DEFAULT_PRINT_FORMAT); };
    // void printf(const char *prefix, const char *format) { printf(prefix, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_COUNT, format); };
    // void printf(const char *prefix) { printf(prefix, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_FORMAT); };
    // void printf(INTE_TYPE row_limit, INTE_TYPE col_limit) { printf("", row_limit, col_limit, COLMAT_DEFAULT_PRINT_FORMAT); };
    // void printf() { printf("", COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_FORMAT); };

    void Ones()
    {
        if (ld == r)
            for (auto ind = 0; ind < r * c; ind++)
                v[ind] = 1;
        else
        {
            ELEM_TYPE *col_ptr = v;
            for (auto col_ind = 0; col_ind < c; col_ind++, col_ptr += ld)
                for (auto row_ind = 0; row_ind < r; row_ind++)
                    col_ptr[row_ind] = 1;
        }
    };

    void Zero()
    {
        if (ld == r)
            memset(v, 0, sizeof(ELEM_TYPE) * r * c);
        else
        {
            ELEM_TYPE *col_ptr = v;
            for (auto col_ind = 0; col_ind < c; col_ind++, col_ptr += ld)
                memset(col_ptr, 0, sizeof(ELEM_TYPE) * r);
        }
    };

    void Scale(REAL_TYPE s)
    {
        if (s == 0)
            Zero();
        else if (ld == r)
            for (auto ind = 0; ind < r * c; ind++)
                v[ind] *= s;
        else
            for (auto col_ind = 0; col_ind < c; col_ind++)
                for (auto row_ind = 0; row_ind < r; row_ind++)
                    v[row_ind + col_ind * ld] *= s;
    }

    void Idmt()
    {
        this->Zero();
        INTE_TYPE minrc = r < c ? r : c;
        for (auto ind = 0; ind < minrc; ind++)
            v[ind + ind * ld] = 1;
    };

    template <class random_engine>
    void Rand(random_engine &engine, ELEM_TYPE a, ELEM_TYPE b)
    {
        std::uniform_real_distribution<ELEM_TYPE> distribution(a, b);

        if (ld == r)
            for (auto ind = 0; ind < r * c; ind++)
                v[ind] = distribution(engine);
        else
        {
            ELEM_TYPE *col_ptr = v;
            for (auto col_ind = 0; col_ind < c; col_ind++, col_ptr += ld)
                for (auto row_ind = 0; row_ind < r; row_ind++)
                    col_ptr[row_ind] = distribution(engine);
        }
    }

    REAL_TYPE Norm1()
    {
        REAL_TYPE norm = 0.0;
        REAL_TYPE sum = 0.0;
        for (auto col_ind = 0; col_ind < c; col_ind++)
        {
            sum = cblas_dasum(r, v + col_ind * ld, 1);
            if (norm < sum)
                norm = sum;
        }
        return norm;
    };
};

template <class ELEM_TYPE>
class ColMat
{
    // The View object NEVER owns the array it accesses. Any View object must be initialized with the array it accesses.
    // It is encouraged to create and release View objects on demand for one-time-use, i.e., it is not recommended to reuse it,
    // i.e., the copy constructor and the assignment constructor are not recommended, as they are called only in reuse scenarios.
    // Also note that ColMat that does not own the vector v, i.e., b = false, functions exactly as View_ColMat with ld = r.

    typedef ColMat<ELEM_TYPE> SELF_TYPE;
    typedef View_ColMat<ELEM_TYPE> VIEW_TYPE;

public:
    ELEM_TYPE *v;
    INTE_TYPE r;
    INTE_TYPE c;

    ColMat() : v(nullptr), r(0), c(0) {};
    ColMat(INTE_TYPE row, INTE_TYPE col) : v(new ELEM_TYPE[row * col]()), r(row), c(col) {};
    void copy(const SELF_TYPE &src)
    {
        // assert(this->r == src.r);
        // assert(this->c == src.c);
        memcpy(this->v, src.v, sizeof(ELEM_TYPE) * r * c);
    };
    ColMat(const SELF_TYPE &src) : ColMat(src.r, src.c) { this->copy(src); };
    void swap(SELF_TYPE &rhs) noexcept
    {
        using std::swap;
        swap(this->v, rhs.v);
        swap(this->r, rhs.r);
        swap(this->c, rhs.c);
    }
    ColMat(SELF_TYPE &&rhs) noexcept : v(rhs.v), r(rhs.r), c(rhs.c)
    {
        rhs.r = 0;
        rhs.c = 0;
        rhs.v = nullptr;
    }
    ColMat<ELEM_TYPE> &operator=(SELF_TYPE rhs) noexcept
    {
        swap(rhs);
        return (*this);
    }
    ~ColMat() { delete[] v; };

    void Assign(ELEM_TYPE *src, INTE_TYPE lds, INTE_TYPE rs, INTE_TYPE cs)
    {

        if ((rs == lds) && (r == lds))
            memcpy(v, src, sizeof(ELEM_TYPE) * rs * cs);
        else
        {
            ELEM_TYPE *des_ptr = v;
            ELEM_TYPE *src_ptr = src;
            for (auto col_ind = 0; col_ind < cs; col_ind++, des_ptr += r, src_ptr += lds)
                memcpy(des_ptr, src_ptr, sizeof(ELEM_TYPE) * rs);
        }
    }
    void Assign(const SELF_TYPE &src) { Assign(src.v, src.r, src.r, src.c); };
    void Assign(const VIEW_TYPE &src) { Assign(src.v, src.ld, src.r, src.c); };
    void Assign(ELEM_TYPE *src, INTE_TYPE lds) { Assign(src, lds, r, c); };
    void Assign(ELEM_TYPE *src) { Assign(src, r, r, c); };

    void Copyto(ELEM_TYPE *des, INTE_TYPE ldd, INTE_TYPE rd, INTE_TYPE cd)
    {
        if ((rd == ldd) && (r == ldd))
            memcpy(des, this->v, sizeof(ELEM_TYPE) * rd * cd);
        else
        {
            ELEM_TYPE *des_col = des;
            ELEM_TYPE *src_col = this->v;
            for (auto col_ind = 0; col_ind < cd; col_ind++, src_col += r, des_col += ldd)
                memcpy(des_col, src_col, sizeof(ELEM_TYPE) * rd);
        }
    };
    void Copyto(const SELF_TYPE &des) { Copyto(des.v, des.r, des.r, des.c); };
    void Copyto(const VIEW_TYPE &des) { Assign(des.v, des.ld, des.r, des.c); };
    void Copyto(ELEM_TYPE *des, INTE_TYPE ldd) { Copyto(des, ldd, r, c); };
    void Copyto(ELEM_TYPE *des) { Copyto(des, r, r, c); };

    ELEM_TYPE &operator()(INTE_TYPE i, INTE_TYPE j) { return *(v + i + j * r); }; // return writable i,j entry
    ELEM_TYPE *operator[](INTE_TYPE i) { return v + i * r; };                     // return the column vector

    ELEM_TYPE *ptr(INTE_TYPE i) { return v + i * r; };
    ELEM_TYPE *ptr() { return v; };

    ELEM_TYPE &cstele(INTE_TYPE i, INTE_TYPE j) const { return *(v + i + j * r); };
    ELEM_TYPE *cstptr(INTE_TYPE i, INTE_TYPE j) const { return v + i + j * r; };

    // It is very important to note that the (deep) copy operator of the View object is NOT called in the copy constructor
    // as this object MUST be an shallow object that does not own or operate any data unless specifically instructed.

    // void fprintf(FILE *of, const char *prefix, INTE_TYPE row_limit, INTE_TYPE col_limit, const char *format)
    // {
    //     using std::fprintf;

    //     INTE_TYPE row = (row_limit < r) ? row_limit : r;
    //     INTE_TYPE col = (col_limit < c) ? col_limit : c;

    //     fprintf(of, "%s", prefix);
    //     fprintf(of, "[\n");
    //     for (auto row_ind = 0; row_ind < row; row_ind++, (col == c ? fprintf(of, "%s", ";\n") : fprintf(of, "%s", "...\n")))
    //         for (auto col_ind = 0; col_ind < col; col_ind++)
    //             fprintf(of, format, ptr(col_ind)[row_ind]);
    //     if (row != r)
    //     {
    //         for (auto col_ind = 0; col_ind < col; col_ind++)
    //             fprintf(of, "%s", "...\t");
    //         fprintf(of, "%s", "\n");
    //     }
    //     fprintf(of, "%s", "]\n\n");
    // };

    // void fprintf(FILE *of, const char *prefix, INTE_TYPE row_limit, INTE_TYPE col_limit) { fprintf(of, prefix, row_limit, col_limit, COLMAT_DEFAULT_PRINT_FORMAT); };
    // void fprintf(FILE *of, const char *prefix, const char *format) { fprintf(of, prefix, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_COUNT, format); };
    // void fprintf(FILE *of, const char *prefix) { fprintf(of, prefix, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_FORMAT); };
    // void fprintf(FILE *of, INTE_TYPE row_limit, INTE_TYPE col_limit) { fprintf(of, "", row_limit, col_limit, COLMAT_DEFAULT_PRINT_FORMAT); };
    // void fprintf(FILE *of) { fprintf(of, "", COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_FORMAT); };

    // void printf(const char *prefix, INTE_TYPE row_limit, INTE_TYPE col_limit, const char *format)
    // {
    //     using std::printf;

    //     INTE_TYPE row = (row_limit < r) ? row_limit : r;
    //     INTE_TYPE col = (col_limit < c) ? col_limit : c;

    //     printf("%s", prefix);
    //     printf("[\n");
    //     for (auto row_ind = 0; row_ind < row; row_ind++, ((col == c) ? printf(";\n") : printf("...\n")))
    //         for (auto col_ind = 0; col_ind < col; col_ind++)
    //             printf(format, ptr(col_ind)[row_ind]);
    //     if (row != r)
    //     {
    //         for (auto col_ind = 0; col_ind < col; col_ind++)
    //             printf("...\t");
    //         printf("\n");
    //     }
    //     printf("]\n\n");
    // };

    // void printf(const char *prefix, INTE_TYPE row_limit, INTE_TYPE col_limit) { printf(prefix, row_limit, col_limit, COLMAT_DEFAULT_PRINT_FORMAT); };
    // void printf(const char *prefix, const char *format) { printf(prefix, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_COUNT, format); };
    // void printf(const char *prefix) { printf(prefix, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_FORMAT); };
    // void printf(INTE_TYPE row_limit, INTE_TYPE col_limit) { printf("", row_limit, col_limit, COLMAT_DEFAULT_PRINT_FORMAT); };
    // void printf() { printf("", COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_COUNT, COLMAT_DEFAULT_PRINT_FORMAT); };

    void Ones()
    {
        for (auto ind = 0; ind < r * c; ind++)
            v[ind] = 1;
    };

    void Zero() { memset(v, 0, sizeof(ELEM_TYPE) * r * c); };

    void Scale(REAL_TYPE s)
    {
        if (s == 0)
            Zero();
        else
            for (auto ind = 0; ind < r * c; ind++)
                v[ind] *= s;
    }

    void Idmt()
    {
        this->Zero();
        INTE_TYPE minrc = r < c ? r : c;
        for (auto ind = 0; ind < minrc; ind++)
            v[ind + ind * r] = 1;
    };

    template <class random_engine>
    void Rand(random_engine &engine, ELEM_TYPE a, ELEM_TYPE b)
    {
        std::uniform_real_distribution<ELEM_TYPE> distribution(a, b);
        for (auto ind = 0; ind < r * c; ind++)
            v[ind] = distribution(engine);
    }

    REAL_TYPE Norm1()
    {
        REAL_TYPE norm = 0.0;
        REAL_TYPE sum = 0.0;
        for (auto col_ind = 0; col_ind < c; col_ind++)
        {
            sum = cblas_dasum(r, v + col_ind * r, 1);
            if (norm < sum)
                norm = sum;
        }
        return norm;
    };
};
