#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <algorithm>
#include "type_convention.hpp"

// The ColMat implemented in this header represents a dense matrix,
// stored in a continuous memory in thecolumn major order.
// The ColMat may or may not own the array of the entries in the matrix.
// The View_ColMat is designed for accessing the sub-matrix in ColMat,
// including itself. Therefore, implicit cast from ColMat to View_ColMat
// is also implemented.

// Only simple data primitives like INTE_TYPE, CMPX_TYPE are supported.

#define MATRIX_INFINITE_NORM -1

template <class type>
class View_ColMat
{
    // The View object NEVER owns the array it is accessing. Any View object must be initialized with the array it is accessing.
    // It is encouraged to create and release View objects on demand for one-time-use, i.e., it is not recommended to reuse it,
    // i.e., the copy constructor and the assignment constructor are not recommended, as they are called only in reuse scenarios.
    // Also note that ColMat that does not own the vector v, i.e., b = false, functions exactly as View_ColMat with ld = r.
public:
    type *v;
    INTE_TYPE r;
    INTE_TYPE c;
    BOOL_TYPE ld;       // The leading order of the view, which is the row number of the underlying ColMat or the equivalent structure.
    type **_col_access; // Fast access to the columns

    View_ColMat() : v(nullptr), r(0), c(0), ld(0), _col_access(nullptr){};
    // View_ColMat(const ColMat<type>& src, INTE_TYPE row, INTE_TYPE col) : v(src.fpele(0,0)), r(row), c(col), ld(src.row()), _col_access(nullptr){};
    // View_ColMat(const ColMat<type>& src) : View_ColMat(src, src.row(), src.col()){};
    // These constructors are replaced by the conversion in ColMat
    View_ColMat(type *vec, INTE_TYPE leading_order, INTE_TYPE row, INTE_TYPE col) : v(vec), r(row), c(col), ld(leading_order), _col_access(nullptr){};
    View_ColMat(const View_ColMat<type> &src) : v(src.v), r(src.r), c(src.c), ld(src.ld)
    {
        if (src._col_access)
        {
            this->_col_access = new type *[c];
            for (INTE_TYPE i = 0; i < c; i++)
                this->_col_access[i] = src._col_access[i];
        }
        else
            this->_col_access = nullptr;
    }

    void swap(View_ColMat<type> &rhs)
    {
        using std::swap;
        swap(this->v, rhs.v);
        swap(this->r, rhs.r);
        swap(this->c, rhs.c);
        swap(this->ld, rhs.ld);
        swap(this->_col_access, rhs._col_access);
    }
    View_ColMat<type> &operator=(const View_ColMat<type> &rhs)
    {
        View_ColMat<type> temp(rhs);
        swap(temp);
        return (*this);
    }
    ~View_ColMat()
    {
        if (_col_access)
            delete[] _col_access;
        _col_access = nullptr;
    }

    void fast_col_access()
    {
        if (!_col_access)
        {
            _col_access = new type *[c];
            for (INTE_TYPE i = 0; i < c; i++)
                _col_access[i] = v + i * ld;
        }
    }

    void release_col_access()
    {
        if (_col_access)
            delete[] _col_access;
        _col_access = nullptr;
    };

    INTE_TYPE nrow() const { return r; };
    INTE_TYPE ncol() const { return c; };

    type ele(INTE_TYPE i, INTE_TYPE j) const { return *(v + i + j * ld); };
    type fele(INTE_TYPE i, INTE_TYPE j) const
    {
        // assert(_col_access != nullptr);
        return *(_col_access[j] + i);
    };
    type *pele(INTE_TYPE i, INTE_TYPE j) const { return v + i + j * ld; };
    type *fpele(INTE_TYPE i, INTE_TYPE j) const
    {
        // assert(_col_access != nullptr);
        return _col_access[j] + i;
    };
    // const access to the entries

    type &operator()(INTE_TYPE i, INTE_TYPE j) { return *(v + i + j * ld); }; // return writable i,j entry
    type *operator[](INTE_TYPE i) { return v + i * ld; };                     // return the column vector
    type &var(INTE_TYPE i, INTE_TYPE j) { return *(v + i + j * ld); };        // same with () operator
    type &fvar(INTE_TYPE i, INTE_TYPE j)
    { // same with () operator, fast access with _col_access
        // assert(_col_access != nullptr);
        return *(_col_access[j] + i);
    };
    type *pvar(INTE_TYPE i, INTE_TYPE j) { return v + i + j * ld; };
    type *fpvar(INTE_TYPE i, INTE_TYPE j)
    {
        // assert(_col_access != nullptr);
        return _col_access[j] + i;
    };
    // writable access to the entries

    type *col(INTE_TYPE i) { return v + i * ld; }
    type *fcol(INTE_TYPE i)
    {
        // assert(_col_access != nullptr);
        return _col_access[i];
    };
    // writable access to column vector

    void reset_array(type *vec, INTE_TYPE leading_order, INTE_TYPE row, INTE_TYPE col)
    {
        // This routine basically functions as the re-assignment to lhs = View_ColMat<type>(vec, leading_order, row, col)
        // but it avoids the multiple constructor, thanks to the simple data structure.
        release_col_access();
        v = vec;
        ld = leading_order;
        r = row;
        c = col;
    };

    void copy(const View_ColMat<type> &src)
    {
        if ((this->r == this->ld) && (src.r == src.ld))
            memcpy(this->v, src.v, sizeof(type) * r * c);
        else
        {
            auto *des_ptr = v;
            auto *src_ptr = src.v;
            for (INTE_TYPE col_ind = 0; col_ind < c; col_ind++, des_ptr += this->ld, src_ptr += src.ld)
                memcpy(des_ptr, src_ptr, sizeof(type) * r);
        }
    }
    // It is very important to note that the (deep) copy operator of the View object is NOT called in the copy constructor
    // as this object MUST be an shallow object that does not own or operate any data unless specifically instructed.

    void assign(const type *src, INTE_TYPE lds)
    {
        type *des_ptr = v;
        type *src_ptr = src;
        for (INTE_TYPE col_ind = 0; col_ind < c; col_ind++, des_ptr += this->ld, src_ptr += lds)
            memcpy(des_ptr, src_ptr, sizeof(type) * r);
    }
    void assign(const type *src)
    {
        if (r == ld)
            memcpy(this->v, src, sizeof(type) * r * c);
        else
        {
            auto *des_ptr = v;
            auto *src_ptr = src;
            for (INTE_TYPE col_ind = 0; col_ind < c; col_ind++, des_ptr += this->ld, src_ptr += this->r)
                memcpy(des_ptr, src_ptr, sizeof(type) * r);
        }
    }
    void assign(const View_ColMat<type> &src)
    {
        this->copy(src);
    }

    void fprintf(FILE *of, INTE_TYPE row_limit, INTE_TYPE col_limit)
    {
        using std::fprintf;

        INTE_TYPE row = (row_limit < r) ? row_limit : r;
        INTE_TYPE col = (col_limit < c) ? col_limit : c;

        for (INTE_TYPE row_ind = 0; row_ind < row; row_ind++, (row_limit < r ? fprintf(of, "...\n") : fprintf(of, "\n")))
            for (INTE_TYPE col_ind = 0; col_ind < col; col_ind++)
                fprintf(of, "%1.3f\t", ele(row_ind, col_ind));
        if (col_limit < c)
        {
            for (INTE_TYPE col_ind = 0; col_ind < col + 1; col_ind++)
                fprintf(of, "...\t");
            fprintf(of, "\n");
        }
        fprintf(of, "\n");
    };
    void printf(INTE_TYPE row_limit, INTE_TYPE col_limit)
    {
        using std::printf;
        INTE_TYPE row = (row_limit < r) ? row_limit : r;
        INTE_TYPE col = (col_limit < c) ? col_limit : c;

        for (INTE_TYPE row_ind = 0; row_ind < row; row_ind++, (row_limit < r ? printf("...\n") : printf("\n")))
            for (INTE_TYPE col_ind = 0; col_ind < col; col_ind++)
                printf("%1.3f\t", ele(row_ind, col_ind));
        if (col_limit < c)
        {
            for (INTE_TYPE col_ind = 0; col_ind < col + 1; col_ind++)
                printf("...\t");
            printf("\n");
        }
        printf("\n");
    };
    void printf() { printf(10, 10); };
};

template <class type>
class ColMat
{
public:
    type *v;
    INTE_TYPE r;
    INTE_TYPE c;
    BOOL_TYPE b; // indicate whether the array v beldng to the ColMat
    type **_col_access;

    // The _col_access should be an on-demand object that needs to be declared whenever needed.
    // Deep copy/Assignment operator/Copy constructor would not initial this object nor copy pointer.
    // Although it is possible to copying _col_access by checking if the source has it initialized,
    // there are concerns on its perfomance impacts.

    ColMat() : v(nullptr), r(0), c(0), b(false), _col_access(nullptr){};
    ColMat(INTE_TYPE row, INTE_TYPE col) : v(new type[row * col]), r(row), c(col), b(true), _col_access(nullptr){};
    ColMat(type *vec, INTE_TYPE row, INTE_TYPE col) : v(vec), r(row), c(col), b(false), _col_access(nullptr){};
    ColMat(const ColMat<type> &src) : ColMat<type>(src.r, src.c) { memcpy(this->v, src.v, sizeof(type) * r * c); };
    void swap(ColMat<type> &rhs)
    {
        using std::swap;
        swap(this->v, rhs.v);
        swap(this->r, rhs.r);
        swap(this->c, rhs.c);
        swap(this->b, rhs.b);
        swap(this->_col_access, rhs._col_access);
    }
    ColMat<type> &operator=(const ColMat<type> &rhs)
    {
        ColMat<type> temp(rhs);
        swap(temp);
        return (*this);
    }
    ~ColMat()
    {
        if (v && b)
            delete[] v;
        v = nullptr;
        if (_col_access)
            delete[] _col_access;
        _col_access = nullptr;
    }
    // For ColMat constructed with v passed in argument, the array v does not beldng to ColMat
    // and it won't be released while ColMat is being released. User or the other struct is responible to manage this array.

    // operator View_ColMat<type>&()
    // {
    //     View_ColMat<type> view(this->v, this->r, this->r, this->c);
    //     return view;
    // }

    operator View_ColMat<type>()
    {
        View_ColMat<type> view(this->v, this->r, this->r, this->c);
        return view;
    };

    //  The implict casting from ColMat to View_ColMat

    void copy(const ColMat<type> &src)
    {
        memcpy(this->v, src.v, sizeof(type) * src.r * src.c);
    };

    void copy(const View_ColMat<type> &src)
    {
        if (src.ld == src.r)
            memcpy(this->v, src.v, sizeof(type) * r * c);
        else
        {
            auto des_ptr = this->v;
            auto src_ptr = src.v;

            for (INTE_TYPE col_ind = 0; col_ind < c; col_ind++, des_ptr += r, src_ptr += src.ld)
                memcpy(des_ptr, src_ptr, sizeof(type) * r);
        }
    };

    void release_col_access()
    {
        if (_col_access)
            delete[] _col_access;
        _col_access = nullptr;
    }

    void fast_col_access()
    {
        if (!_col_access)
        {
            _col_access = new type *[c];
            for (INTE_TYPE i = 0; i < c; i++)
                _col_access[i] = v + i * r;
        }
    }

    INTE_TYPE nrow() const { return r; };
    INTE_TYPE ncol() const { return c; };

    type ele(INTE_TYPE i, INTE_TYPE j) const { return *(v + i + j * r); };
    type fele(INTE_TYPE i, INTE_TYPE j) const
    {
        // assert(_col_access != nullptr);
        return *(_col_access[j] + i);
    };
    type *pele(INTE_TYPE i, INTE_TYPE j) const { return v + i + j * r; };
    type *fpele(INTE_TYPE i, INTE_TYPE j) const
    {
        // assert(_col_access != nullptr);
        return _col_access[j] + i;
    };
    // const access to the entries

    type &operator()(INTE_TYPE i, INTE_TYPE j) { return *(v + i + j * r); }; // return writable i,j entry
    type *operator[](INTE_TYPE i) { return v + i * r; };                     // return the column vector
    type &var(INTE_TYPE i, INTE_TYPE j) { return *(v + i + j * r); };        // same with () operator
    type &fvar(INTE_TYPE i, INTE_TYPE j)
    { // same with () operator, fast access with _col_access
        // assert(_col_access != nullptr);
        return *(_col_access[j] + i);
    };
    type *pvar(INTE_TYPE i, INTE_TYPE j) { return v + i + j * r; };
    type *fpvar(INTE_TYPE i, INTE_TYPE j)
    {
        // assert(_col_access != nullptr);
        return _col_access[j] + i;
    };
    // writable access to the entries

    type *col(INTE_TYPE i) { return v + i * r; }
    type *fcol(INTE_TYPE i)
    {
        assert(_col_access != nullptr);
        return _col_access[i];
    };
    // writable access to column vector

    BOOL_TYPE vec_owner() const { return b; };

    void reset_array(type *vec, INTE_TYPE row, INTE_TYPE col)
    {
        // This routine basically functions as the re-assignment to lhs = ColMat<type>(vec, row, col)
        // but it avoids the multiple constructor, thanks to the simple data structure.
        // assert(b != true);
        // This function must never be used for ColMat that is responible to the array v.
        release_col_access();
        b = false;
        v = vec;
        r = row;
        c = col;
    }

    REAL_TYPE opnorm(INTE_TYPE p)
    {
        REAL_TYPE norm = 0.0;
        REAL_TYPE col_norm = 0.0;
        type *ptr = this->v;
        if (p == 1)
        {
            for (INTE_TYPE col_ind = 0; col_ind < this->c; col_ind++, ptr += this->r, norm = (col_norm > norm) ? col_norm : norm, col_norm = 0.0)
                for (INTE_TYPE row_ind = 0; row_ind < this->r; row_ind++)
                    col_norm += abs(*(ptr + row_ind));
        }
        else
        {
            std::printf("Vector norm %lld not supported!", p);
            throw(1);
        }

        return norm;
    }

    void Skew()
    {
        // assert(r == c);
        REAL_TYPE temp = 0.0;
        if (_col_access)
            for (INTE_TYPE col_ind = 0; col_ind < c; col_ind++)
                for (INTE_TYPE row_ind = col_ind; row_ind < r; fvar(col_ind, row_ind) = -fvar(row_ind, col_ind), row_ind++)
                    fvar(row_ind, col_ind) = (fvar(row_ind, col_ind) - fvar(col_ind, row_ind)) * 0.5;
        else
            for (INTE_TYPE col_ind = 0; col_ind < c; col_ind++)
                for (INTE_TYPE row_ind = col_ind; row_ind < r; v[col_ind + row_ind * r] = -v[row_ind + col_ind * r], row_ind++)
                    v[row_ind + col_ind * r] = (v[row_ind + col_ind * r] - v[col_ind + row_ind * r]) * 0.5;
    };

    void Symm()
    {
        // assert(r == c);
        REAL_TYPE temp = 0.0;
        if (_col_access)
            for (INTE_TYPE col_ind = 0; col_ind < c; col_ind++)
                for (INTE_TYPE row_ind = col_ind; row_ind < r; fvar(col_ind, row_ind) = fvar(row_ind, col_ind), row_ind++)
                    fvar(row_ind, col_ind) = (fvar(row_ind, col_ind) + fvar(col_ind, row_ind)) * 0.5;
        else
            for (INTE_TYPE col_ind = 0; col_ind < c; col_ind++)
                for (INTE_TYPE row_ind = col_ind; row_ind < r; v[col_ind + row_ind * r] = v[row_ind + col_ind * r], row_ind++)
                    v[row_ind + col_ind * r] = (v[row_ind + col_ind * r] + v[col_ind + row_ind * r]) * 0.5;
    };

    void fprintf(FILE *of, INTE_TYPE row_limit, INTE_TYPE col_limit)
    {
        using std::fprintf;

        INTE_TYPE row = (row_limit < r) ? row_limit : r;
        INTE_TYPE col = (col_limit < c) ? col_limit : c;

        for (INTE_TYPE row_ind = 0; row_ind < row; row_ind++, (row_limit < r ? fprintf(of, "...\n") : fprintf(of, "\n")))
            for (INTE_TYPE col_ind = 0; col_ind < col; col_ind++)
                fprintf(of, "%1.3f\t", ele(row_ind, col_ind));
        if (col_limit < c)
        {
            for (INTE_TYPE col_ind = 0; col_ind < col + 1; col_ind++)
                fprintf(of, "...\t");
            fprintf(of, "\n");
        }
        fprintf(of, "\n");
    };
    void printf(INTE_TYPE row_limit, INTE_TYPE col_limit)
    {
        using std::printf;
        INTE_TYPE row = (row_limit < r) ? row_limit : r;
        INTE_TYPE col = (col_limit < c) ? col_limit : c;

        for (INTE_TYPE row_ind = 0; row_ind < row; row_ind++, (row_limit < r ? printf("...\n") : printf("\n")))
            for (INTE_TYPE col_ind = 0; col_ind < col; col_ind++)
                printf("%1.3f\t", ele(row_ind, col_ind));
        if (col_limit < c)
        {
            for (INTE_TYPE col_ind = 0; col_ind < col + 1; col_ind++)
                printf("...\t");
            printf("\n");
        }
        printf("\n");
    };
    void printf() { printf(10, 10); };
    void printf(const char *s, INTE_TYPE row_limit, INTE_TYPE col_limit)
    {
        std::cout << s;
        printf(row_limit, col_limit);
    };
    void printf(const char *s)
    {
        std::cout << s;
        printf(10, 10);
    };

    void assign(const type *src, INTE_TYPE lds)
    {
        auto *des_ptr = v;
        auto *src_ptr = src;
        for (INTE_TYPE col_ind = 0; col_ind < c; col_ind++, des_ptr += r, src_ptr += lds)
            memcpy(des_ptr, src_ptr, sizeof(type) * r);
    }
    void assign(const type *src)
    {
        memcpy(v, src, sizeof(type) * r * c);
    }
    void assign(const ColMat<type> &src)
    {
        this->copy(src);
    }
    void assign(const View_ColMat<type> &src)
    {
        this->copy(src);
    }

    // friend void copy(const View_ColMat<type> &ViewDes, const ColMat<type> &MatSrc)
    //{
    //     if (ViewDes._col_access)
    //         if (MatSrc._col_access)
    //             for (INTE_TYPE col_ind = 0; col_ind < MatSrc.c; col_ind++)
    //                 memcpy(ViewDes._col_access[col_ind], MatSrc._col_access[col_ind], sizeof(type) * MatSrc.r);
    //         else
    //             for (INTE_TYPE col_ind = 0; col_ind < MatSrc.c; col_ind++)
    //                 memcpy(ViewDes._col_access[col_ind], MatSrc.v + (col_ind * MatSrc.r), sizeof(type) * MatSrc.r);
    //     else if (MatSrc._col_access)
    //         for (INTE_TYPE col_ind = 0; col_ind < MatSrc.c; col_ind++)
    //             memcpy(ViewDes.v + (col_ind * ViewDes.ld), MatSrc._col_access[col_ind], sizeof(type) * MatSrc.r);
    //     else
    //         for (INTE_TYPE col_ind = 0; col_ind < MatSrc.c; col_ind++)
    //             memcpy(ViewDes.v + (col_ind * ViewDes.ld), MatSrc.v + (col_ind * MatSrc.r), sizeof(type) * MatSrc.r);
    // };
};
