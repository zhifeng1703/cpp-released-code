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
// #include "mkl.h"

#define ARRVEC_DEFAULT_PRINT_FORMAT "%1.3f\t"
#ifndef ARRVEC_DEFAULT_PRINT_COUNT
#define ARRVEC_DEFAULT_PRINT_COUNT 10
#endif

template <class ELEM_TYPE>
class View_ArrVec
{
    typedef View_ArrVec<ELEM_TYPE> SELF_TYPE;

public:
    ELEM_TYPE *v;
    INTE_TYPE d;

    View_ArrVec() : v(nullptr), d(0) {};
    View_ArrVec(ELEM_TYPE *vec, INTE_TYPE dim) : v(vec), d(dim) {};
    void copy(const SELF_TYPE &src)
    {
        // no deep copy should be executed for View ELEM_TYPE.
        this->v = src.v;
        this->d = src.d;
    };
    View_ArrVec(const SELF_TYPE &src) : View_ArrVec<ELEM_TYPE>(src.v, src.d) {};
    void swap(SELF_TYPE &rhs)
    {
        std::swap(this->v, rhs.v);
        std::swap(this->d, rhs.d);
    }
    SELF_TYPE &operator=(const SELF_TYPE &rhs)
    {
        SELF_TYPE temp(rhs);
        swap(temp);
        return (*this);
    }
    ~View_ArrVec() {}; // View ELEM_TYPE should never release any memory it accesses.

    operator ELEM_TYPE *() { return v; };

    void Assign(ELEM_TYPE *src, INTE_TYPE length) { memcpy(this->v, src, sizeof(ELEM_TYPE) * length); };
    void Assign(const SELF_TYPE &src) { Assign(src.v, src.d); };
    void Assign(ELEM_TYPE *src) { Assign(src, this->d); };

    void Copyto(ELEM_TYPE *des, INTE_TYPE length) { memcpy(des, this->v, sizeof(ELEM_TYPE) * length); };
    void Copyto(const SELF_TYPE &des) { Copyto(des.v, this->d); };
    void Copyto(ELEM_TYPE *des) { Copyto(des, this->d); };

    ELEM_TYPE &operator()(INTE_TYPE i) { return *(v + i); };
    ELEM_TYPE &operator[](INTE_TYPE i) { return *(v + i); };

    ELEM_TYPE *ptr(INTE_TYPE i) { return v + i; };
    ELEM_TYPE *ptr() { return v; };

    void Ones()
    {
        for (auto vec_ind = 0; vec_ind < d; vec_ind++)
            v[vec_ind] = 1;
    };
    void Zero() { memset(v, 0, sizeof(ELEM_TYPE) * d); };

    template <class random_engine>
    void Rand(random_engine &engine, ELEM_TYPE a, ELEM_TYPE b)
    {
        std::uniform_real_distribution<ELEM_TYPE> distribution(a, b);
        for (auto vec_ind = 0; vec_ind < d; vec_ind++)
            v[vec_ind] = distribution(engine);
    }

    ELEM_TYPE &cstele(INTE_TYPE i) const { return *(v + i); };
    ELEM_TYPE *cstptr(INTE_TYPE i) const { return v + i; };
    // const access to the entries

    void Reset(ELEM_TYPE *vec, INTE_TYPE dim)
    {
        v = vec;
        d = dim;
    };

    // void fprintf(FILE *of, const char *prefix, INTE_TYPE length, const char *format)
    // {
    //     using std::fprintf;
    //     INTE_TYPE dim = (length < d) ? length : d;

    //     fprintf(of, prefix);
    //     for (auto vec_ind = 0; vec_ind < dim; vec_ind++)
    //         fprintf(of, format, v[vec_ind]);
    //     if (length < d)
    //         fprintf(of, "...\n\n");
    //     else
    //         fprintf(of, "\n\n");
    // };

    // void fprintf(FILE *of, const char *prefix, INTE_TYPE length) { fprintf(of, prefix, length, ARRVEC_DEFAULT_PRINT_FORMAT); };
    // void fprintf(FILE *of, const char *prefix, const char *format) { fprintf(of, prefix, ARRVEC_DEFAULT_PRINT_COUNT, format); };
    // void fprintf(FILE *of, const char *prefix) { fprintf(of, prefix, ARRVEC_DEFAULT_PRINT_COUNT, ARRVEC_DEFAULT_PRINT_FORMAT); };

    // void fprintf(FILE *of, INTE_TYPE length, const char *format) { fprintf(of, "", length, format); };
    // void fprintf(FILE *of, INTE_TYPE length) { fprintf(of, length, ARRVEC_DEFAULT_PRINT_FORMAT); };
    // void fprintf(FILE *of) { fprintf(of, ARRVEC_DEFAULT_PRINT_COUNT, ARRVEC_DEFAULT_PRINT_FORMAT); };

    // void printf(const char *prefix, INTE_TYPE length, const char *format)
    // {
    //     using std::printf;
    //     INTE_TYPE dim = (length < d) ? length : d;
    //     printf("%s", prefix);
    //     for (auto vec_ind = 0; vec_ind < dim; vec_ind++)
    //         printf(format, v[vec_ind]);
    //     if (length < d)
    //         printf("...\n\n");
    //     else
    //         printf("\n\n");
    // };

    // void printf(const char *prefix, INTE_TYPE length) { printf(prefix, length, ARRVEC_DEFAULT_PRINT_FORMAT); };
    // void printf(const char *prefix, const char *format) { printf(prefix, ARRVEC_DEFAULT_PRINT_COUNT, format); };
    // void printf(const char *prefix) { printf(prefix, ARRVEC_DEFAULT_PRINT_COUNT, ARRVEC_DEFAULT_PRINT_FORMAT); };

    // void printf(INTE_TYPE length) { printf("", length, ARRVEC_DEFAULT_PRINT_FORMAT); };
    // void printf() { printf("", ARRVEC_DEFAULT_PRINT_COUNT, ARRVEC_DEFAULT_PRINT_FORMAT); };
};

template <class ELEM_TYPE>
class ArrVec
{
    typedef ArrVec<ELEM_TYPE> SELF_TYPE;

public:
    ELEM_TYPE *v;
    INTE_TYPE d;

    ArrVec() : v(nullptr), d(0) {};
    ArrVec(INTE_TYPE dim) : v(new ELEM_TYPE[dim]()), d(dim) {};
    void copy(const SELF_TYPE &src)
    {
        if (this->d == src.d)
            memcpy(this->v, src.v, sizeof(ELEM_TYPE) * d);
        else
        {
            delete[] v;
            this->d = src.d;
            v = new ELEM_TYPE[d];
            memcpy(this->v, src.v, sizeof(ELEM_TYPE) * d);
        }
    };
    ArrVec(const SELF_TYPE &src) : ArrVec<ELEM_TYPE>(src.d) { this->copy(src); };
    void swap(SELF_TYPE &rhs) noexcept
    {
        using std::swap;
        swap(this->v, rhs.v);
        swap(this->d, rhs.d);
    }
    ArrVec(SELF_TYPE &&rhs) noexcept : v(rhs.v), d(rhs.d)
    {
        rhs.d = 0;
        rhs.v = nullptr;
    }
    SELF_TYPE &operator=(SELF_TYPE rhs) noexcept
    {
        swap(rhs);
        return (*this);
    }
    ~ArrVec() { delete[] v; };

    operator ELEM_TYPE *() const { return v; };
    operator View_ArrVec<ELEM_TYPE>()
    {
        View_ArrVec<ELEM_TYPE> view(this->v, this->d);
        return view;
    };
    //  The casting from ColMat to View_ColMat

    void Assign(ELEM_TYPE *src, INTE_TYPE n) { memcpy(this->v, src, sizeof(ELEM_TYPE) * n); };
    void Assign(ELEM_TYPE *src) { memcpy(this->v, src, sizeof(ELEM_TYPE) * this->d); };
    void Assign(const ArrVec<ELEM_TYPE> &src) { memcpy(this->v, src.v, sizeof(ELEM_TYPE) * src.d); };
    void Assign(const View_ArrVec<ELEM_TYPE> &src) { memcpy(this->v, src.v, sizeof(ELEM_TYPE) * src.d); };

    void Copyto(ELEM_TYPE *des, INTE_TYPE n) { memcpy(des, this->v, sizeof(ELEM_TYPE) * n); };
    void Copyto(ELEM_TYPE *des) { Copyto(des, d); };
    void Copyto(const ArrVec<ELEM_TYPE> &des) { memcpy(des.v, this->v, sizeof(ELEM_TYPE) * this->d); };
    void Copyto(const View_ArrVec<ELEM_TYPE> &des) { memcpy(des.v, this->v, sizeof(ELEM_TYPE) * this->d); };

    ELEM_TYPE &operator()(INTE_TYPE i) { return *(v + i); };
    ELEM_TYPE &operator[](INTE_TYPE i) { return *(v + i); };

    ELEM_TYPE *ptr(INTE_TYPE i) { return v + i; };
    ELEM_TYPE *ptr() { return v; };

    void Ones(INTE_TYPE n)
    {
        for (auto vec_ind = 0; vec_ind < n; vec_ind++)
            v[vec_ind] = 1;
    };

    void Zero(INTE_TYPE n) { memset(v, 0, sizeof(ELEM_TYPE) * n); };

    template <class random_engine>
    void Rand(INTE_TYPE n, random_engine &engine, ELEM_TYPE a, ELEM_TYPE b)
    {
        std::uniform_real_distribution<ELEM_TYPE> distribution(a, b);
        for (auto vec_ind = 0; vec_ind < n; vec_ind++)
            v[vec_ind] = distribution(engine);
    }

    void Ones() { this->Ones(d); };
    void Zero() { this->Zero(d); };

    template <class random_engine>
    void Rand(random_engine &engine, ELEM_TYPE a, ELEM_TYPE b) { this->Rand(d, engine, a, b); };

    ELEM_TYPE &cstele(INTE_TYPE i) const { return *(v + i); };
    ELEM_TYPE *cstptr(INTE_TYPE i) const { return v + i; };

    // const access to the entries

    // virtual void fprintf(FILE *of, const char *prefix, INTE_TYPE length, const char *format)
    // {
    //     using std::fprintf;
    //     INTE_TYPE dim = (length < d) ? length : d;

    //     fprintf(of, "%s", prefix);
    //     for (auto vec_ind = 0; vec_ind < dim; vec_ind++)
    //         fprintf(of, format, v[vec_ind]);
    //     if (length < d)
    //         fprintf(of, "...\n\n");
    //     else
    //         fprintf(of, "\n\n");
    // };

    // void fprintf(FILE *of, const char *prefix, INTE_TYPE length) { fprintf(of, prefix, length, ARRVEC_DEFAULT_PRINT_FORMAT); };
    // void fprintf(FILE *of, const char *prefix, const char *format) { fprintf(of, prefix, ARRVEC_DEFAULT_PRINT_COUNT, format); };
    // void fprintf(FILE *of, const char *prefix) { fprintf(of, prefix, ARRVEC_DEFAULT_PRINT_COUNT, ARRVEC_DEFAULT_PRINT_FORMAT); };

    // void fprintf(FILE *of, INTE_TYPE length, const char *format) { fprintf(of, "", length, format); };
    // void fprintf(FILE *of, INTE_TYPE length) { fprintf(of, length, ARRVEC_DEFAULT_PRINT_FORMAT); };
    // void fprintf(FILE *of) { fprintf(of, ARRVEC_DEFAULT_PRINT_COUNT, ARRVEC_DEFAULT_PRINT_FORMAT); };

    // virtual void printf(const char *prefix, INTE_TYPE length, const char *format)
    // {
    //     using std::printf;
    //     INTE_TYPE dim = (length < d) ? length : d;
    //     printf("%s", prefix);
    //     for (auto vec_ind = 0; vec_ind < dim; vec_ind++)
    //         printf(format, v[vec_ind]);
    //     if (length < d)
    //         printf("...\n\n");
    //     else
    //         printf("\n\n");
    // };

    // void printf(const char *prefix, INTE_TYPE length) { printf(prefix, length, ARRVEC_DEFAULT_PRINT_FORMAT); };
    // void printf(const char *prefix, const char *format) { printf(prefix, ARRVEC_DEFAULT_PRINT_COUNT, format); };
    // void printf(const char *prefix) { printf(prefix, ARRVEC_DEFAULT_PRINT_COUNT, ARRVEC_DEFAULT_PRINT_FORMAT); };

    // void printf(INTE_TYPE length) { printf("", length, ARRVEC_DEFAULT_PRINT_FORMAT); };
    // void printf() { printf("", ARRVEC_DEFAULT_PRINT_COUNT, ARRVEC_DEFAULT_PRINT_FORMAT); };
};
