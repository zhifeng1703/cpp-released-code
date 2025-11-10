#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>

#include "blasType.hpp"
#include "arrVec.hpp"
#include "colMat.hpp"

// The traversal consist with the following 3 maps. Please follows the guideline to create a customized traversal.
// (1) mat2vec, size: (dim * dim), it maps the columns major matrix indices, (i - 1) + (j - 1) * dim of the (i, j) entry, to the lower-tri vector indices.
//      When (i > j), the entry is on the lower triangular part, return the normal mat2vec((i - 1) + (j - 1) * dim).
//      When (i < j), the entry is on the upper triangular part, return the negation of mat2vec((j - 1) + (i - 1) * dim).
//      When (i == j), the entry on diagonal and sometime does not carry info. Return dim * (dim - 1) / 2 + i at the end of vec that may not exist.
// (2) vec2low, size: (dim * (dim - 1) / 2) + dim, it maps the lower-tri vector indices to the column major matrix indices in the lower triangular part.
// (3) vec2upp, size: (dim * (dim - 1) / 2) + dim, it maps the lower-tri vector indices to the column major matrix indices of the upper triangular part.
// (4) In the case when the diagonal carries no info, the diagonals of mat always map to the last part of the vec, i.e.,
//      mat2vec[i + dim * i] = (dim * (dim - 1) / 2) + i
//      vec2low[(dim * (dim - 1) / 2) + i] = i + dim * i
//      vec2upp[(dim * (dim - 1) / 2) + i] = i + dim * i
// Note that (2) and (3) together essentially store the inverse map of (1).

void strict_lower_blk_traversal(INTE_TYPE *, INTE_TYPE *, INTE_TYPE *, INTE_TYPE);
void strict_lower_col_traversal(INTE_TYPE *, INTE_TYPE *, INTE_TYPE *, INTE_TYPE);

enum LOWTRA_TYPE
{
    BLK_OFFD,
    COL_DIAG,
    COL_OFFD,
    OTHER
};

class LowerTraversal
{
    typedef LowerTraversal SELF_TYPE;

public:
    LOWTRA_TYPE lt_type;
    INTE_TYPE *m2v;
    INTE_TYPE *v2l;
    INTE_TYPE *v2u;
    INTE_TYPE d; // Dimension of the matrix
    INTE_TYPE l; // Size of the vector storing the lower triangular matrix, depending on whether the diagonals are stored.

    LowerTraversal() : lt_type(COL_DIAG), m2v(nullptr), v2l(nullptr), v2u(nullptr), d(0), l(0) {};
    LowerTraversal(INTE_TYPE dim) : lt_type(COL_DIAG), m2v(nullptr), v2l(nullptr), v2u(nullptr), d(dim), l((dim + 1) * dim / 2) {};
    void release()
    {
        if (m2v)
            delete[] m2v;
        if (v2l)
            delete[] v2l;
        if (v2u)
            delete[] v2u;
        m2v = nullptr;
        v2l = nullptr;
        v2u = nullptr;
    }
    void copy(const SELF_TYPE &src)
    {
        lt_type = src.lt_type;
        if (d != src.d)
            this->release();
        d = src.d;
        l = src.l;

        if (src.m2v != nullptr)
        {
            m2v = new INTE_TYPE[d * d];
            memcpy(m2v, src.m2v, sizeof(INTE_TYPE) * d * d);
        }
        if (src.v2l != nullptr)
        {
            v2l = new INTE_TYPE[((d + 1) * d) / 2];      // Space for the diagonal are always assigned but not necessarily used.
            memcpy(v2l, src.v2l, sizeof(INTE_TYPE) * l); // Therefore, the memcopy only copy the actual usage of v2l, in the first l entries.
        }
        if (src.v2u != nullptr)
        {
            v2u = new INTE_TYPE[((d + 1) * d) / 2];      // Space for the diagonal are always assigned but not necessarily used.
            memcpy(v2u, src.v2u, sizeof(INTE_TYPE) * l); // Therefore, the memcopy only copy the actual usage of v2u, in the first l entries.
        }
    }
    LowerTraversal(const SELF_TYPE &src) : LowerTraversal(src.d) { this->copy(src); };
    void swap(SELF_TYPE &rhs)
    {
        using std::swap;
        swap(this->lt_type, rhs.lt_type);
        swap(this->m2v, rhs.m2v);
        swap(this->v2l, rhs.v2l);
        swap(this->v2u, rhs.v2u);
        swap(this->d, rhs.d);
        swap(this->l, rhs.l);
    };
    LowerTraversal(INTE_TYPE dim, LOWTRA_TYPE lt)
    {
        lt_type = lt;
        if (lt == COL_DIAG)
        {
            d = dim;
            l = (d + 1) * d / 2;
            m2v = nullptr;
            v2l = nullptr;
            v2u = nullptr;
        }
        else if (lt == COL_OFFD)
        {
            d = dim;
            l = (d - 1) * d / 2;
            m2v = nullptr;
            v2l = nullptr;
            v2u = nullptr;
        }
        else if (lt == BLK_OFFD)
        {
            d = dim;
            l = (d - 1) * d / 2;
            m2v = new INTE_TYPE[d * d];
            v2l = new INTE_TYPE[((d + 1) * d) / 2];
            v2u = new INTE_TYPE[((d + 1) * d) / 2];
            _blk_offdiagonal_traversal(m2v, v2l, v2u, d);
        }
        else
        {
            std::cout << "Error! LowerTraversal on nonspecified type need to provide construction routine ";
            std::cout << "in the format of void (*traversal)(INTE_TYPE *m2v, INTE_TYPE *v2l, INTE_TYPE *v2u, INTE_TYPE d).\n";
            throw(1);
        }
    }
    LowerTraversal(INTE_TYPE dim, BOOL_TYPE strict_lower, void (*traversal)(INTE_TYPE *, INTE_TYPE *, INTE_TYPE *, INTE_TYPE))
    {
        lt_type = OTHER;
        d = dim;
        l = strict_lower ? ((d - 1) * d) / 2 : ((d + 1) * d) / 2;
        m2v = new INTE_TYPE[d * d];
        v2l = new INTE_TYPE[((d + 1) * d) / 2];
        v2u = new INTE_TYPE[((d + 1) * d) / 2];
        traversal(m2v, v2l, v2u, dim);
    };
    LowerTraversal &operator=(const LowerTraversal &rhs)
    {
        LowerTraversal temp(rhs);
        swap(temp);
        return (*this);
    };
    ~LowerTraversal() { this->release(); };

    template <class ELEM_TYPE>
    void mat2vec(ELEM_TYPE *vec, ELEM_TYPE *mat, INTE_TYPE ld) const
    {
        if (lt_type == OTHER || lt_type == BLK_OFFD) // Nontrivial traversal is called, which is assumed to have the indices maps specified.
            if (ld == d)                             // The leading dimension of mat must be the same with the specified dimension to use v2l array.
                for (auto ind = 0; ind < l; ind++)
                    vec[ind] = mat[v2l[ind]];
            else
            {
                std::cout << "Error! Nontrival traversal by specified indices map is only applicable when the leading dimension of the matrix agrees with the dimension.\n";
                throw(1);
            }
        else if (lt_type == COL_DIAG)
        {
            ELEM_TYPE *col_ptr = mat;
            ELEM_TYPE *vec_ptr = vec;
            for (auto col_ind = 0; col_ind < d; col_ptr += ld + 1, vec_ptr += d - col_ind, col_ind++)
                memcpy(vec_ptr, col_ptr, sizeof(ELEM_TYPE) * (d - col_ind));
        }
        else if (lt_type == COL_OFFD)
        {
            ELEM_TYPE *col_ptr = mat + 1;
            ELEM_TYPE *vec_ptr = vec;
            for (auto col_ind = 0; col_ind < d - 1; col_ptr += ld + 1, vec_ptr += d - col_ind - 1, col_ind++)
                memcpy(vec_ptr, col_ptr, sizeof(ELEM_TYPE) * (d - col_ind - 1));
        }
        else
        {
            std::cout << "Error! Traversal type not supported.\n";
            throw(1);
        }
    }

    template <class ELEM_TYPE>
    void vec2low(ELEM_TYPE *mat, INTE_TYPE ld, ELEM_TYPE *vec) const
    {
        if (lt_type == OTHER || lt_type == BLK_OFFD) // Nontrivial traversal is called, which is assumed to have the indices maps specified.
            if (ld == d)                             // The leading dimension of mat must be the same with the specified dimension to use v2l array.
                for (auto ind = 0; ind < l; ind++)
                    mat[v2l[ind]] = vec[ind];
            else
            {
                std::cout << "Error! Nontrival traversal by specified indices map is only applicable when the leading dimension of the matrix agrees with the dimension.\n";
                throw(1);
            }
        else if (lt_type == COL_DIAG)
        {
            ELEM_TYPE *col_ptr = mat;
            ELEM_TYPE *vec_ptr = vec;
            for (auto col_ind = 0; col_ind < d; col_ptr += ld + 1, vec_ptr += d - col_ind, col_ind++)
                memcpy(col_ptr, vec_ptr, sizeof(ELEM_TYPE) * (d - col_ind));
        }
        else if (lt_type == COL_OFFD)
        {
            ELEM_TYPE *col_ptr = mat + 1;
            ELEM_TYPE *vec_ptr = vec;
            for (auto col_ind = 0; col_ind < d - 1; col_ptr += ld + 1, vec_ptr += d - col_ind - 1, col_ind++)
                memcpy(col_ptr, vec_ptr, sizeof(ELEM_TYPE) * (d - col_ind - 1));
        }
        else
        {
            std::cout << "Error! Traversal type not supported.\n";
            throw(1);
        }
    }

    template <class ELEM_TYPE>
    void vec2upp(ELEM_TYPE *mat, INTE_TYPE ld, ELEM_TYPE *vec) const
    {
        if (lt_type == OTHER || lt_type == BLK_OFFD) // Nontrivial traversal is called, which is assumed to have the indices maps specified.
            if (ld == d)                             // The leading dimension of mat must be the same with the specified dimension to use v2l array.
                for (auto ind = 0; ind < l; ind++)
                    mat[v2u[ind]] = vec[ind];
            else
            {
                std::cout << "Error! Nontrival traversal by specified indices map is only applicable when the leading dimension of the matrix agrees with the dimension.\n";
                throw(1);
            }
        else if (lt_type == COL_DIAG)
        {
            ELEM_TYPE *col_ptr = mat;
            ELEM_TYPE *vec_ptr = vec;
            for (auto col_ind = 0; col_ind < d; col_ptr += ld + 1, vec_ptr += d - col_ind, col_ind++)
                for (auto row_ind = 0; row_ind < d - col_ind; row_ind++)
                    col_ptr[row_ind * ld] = vec_ptr[row_ind];
        }
        else if (lt_type == COL_OFFD)
        {
            ELEM_TYPE *col_ptr = mat + ld;
            ELEM_TYPE *vec_ptr = vec;
            for (auto col_ind = 0; col_ind < d - 1; col_ptr += ld + 1, vec_ptr += d - col_ind - 1, col_ind++)
                for (auto row_ind = 0; row_ind < d - col_ind - 1; row_ind++)
                    col_ptr[row_ind * ld] = vec_ptr[row_ind];
        }
        else
        {
            std::cout << "Error! Traversal type not supported.\n";
            throw(1);
        }
    }

    void _blk_offdiagonal_traversal(INTE_TYPE *mat2vec, INTE_TYPE *vec2low, INTE_TYPE *vec2upp, INTE_TYPE dim)
    {
        INTE_TYPE msize = dim * dim;
        INTE_TYPE lsize = ((dim + 1) * dim) / 2;
        INTE_TYPE bsize = dim / 2;

        INTE_TYPE vec_index = 0;
        for (INTE_TYPE j = 0; j < bsize; j++)
        {
            for (INTE_TYPE i = j + 1; i < bsize; i++)
            {
                // 2 x 2 block at [2i:(2i+1), 2j:(2j+1)], strict lower condition implies i > j.
                // The block is also in column major order, i.e., (1, 0) -> (2, 0) -> (m, 0) -> (2, 1) -> (2, 2) -> ... -> (m, m - 1), where m = bsize - 1
                mat2vec[2 * i + 2 * j * dim] = vec_index;
                vec2low[vec_index] = 2 * i + 2 * j * dim;
                vec2upp[vec_index] = 2 * j + 2 * i * dim;
                vec_index++;

                mat2vec[(2 * i + 1) + 2 * j * dim] = vec_index;
                vec2low[vec_index] = (2 * i + 1) + 2 * j * dim;
                vec2upp[vec_index] = 2 * j + (2 * i + 1) * dim;
                vec_index++;

                mat2vec[2 * i + (2 * j + 1) * dim] = vec_index;
                vec2low[vec_index] = 2 * i + (2 * j + 1) * dim;
                vec2upp[vec_index] = (2 * j + 1) + 2 * i * dim;
                vec_index++;

                mat2vec[(2 * i + 1) + (2 * j + 1) * dim] = vec_index;
                vec2low[vec_index] = (2 * i + 1) + (2 * j + 1) * dim;
                vec2upp[vec_index] = (2 * j + 1) + (2 * i + 1) * dim;
                vec_index++;
            }
        }

        // If there is a left-over row and column.
        if (2 * bsize != dim)
        {
            for (INTE_TYPE j = 0; j < bsize; j++)
            {
                // 1 x 2 block at [dim, 2j:(2j+1)], in the order (dim - 1, 0) -> (dim - 1, 1) -> ... -> (dim - 1, m), where m = bsize - 1.
                mat2vec[(dim - 1) + 2 * j * dim] = vec_index;
                vec2low[vec_index] = (dim - 1) + 2 * j * dim;
                vec2upp[vec_index] = 2 * j + (dim - 1) * dim;
                vec_index++;

                mat2vec[(dim - 1) + (2 * j + 1) * dim] = vec_index;
                vec2low[vec_index] = (dim - 1) + (2 * j + 1) * dim;
                vec2upp[vec_index] = (2 * j + 1) + (dim - 1) * dim;
                vec_index++;
            }
        }

        // Take care of left over entries in the lower left of the diagonal 2 x 2 blocks
        for (INTE_TYPE i = 0; i < bsize; i++)
        {
            // 2 x 2 diagonal blocks in [2i:(2i+1), 2i:(2i+1)] at entry [2i+1, 2i]
            mat2vec[(2 * i + 1) + 2 * i * dim] = vec_index;
            vec2low[vec_index] = (2 * i + 1) + 2 * i * dim;
            vec2upp[vec_index] = 2 * i + (2 * i + 1) * dim;
            vec_index++;
        }

        // Take care of the diagonals
        for (INTE_TYPE i = 0; i < dim; i++)
        {
            mat2vec[i + i * dim] = vec_index;
            vec2low[vec_index] = i + i * dim;
            vec2upp[vec_index] = i + i * dim;
            vec_index++;
        }
    }
};
