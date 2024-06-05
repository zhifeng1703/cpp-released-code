#pragma once

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>

#include "type_convention.hpp"

// The traversal consist with the following 3 maps. Please follows the guideline to create a customized traversal.
// (1) mat2vec, size: (dim * dim), it maps the columns major matrix indices, (i - 1) + (j - 1) * dim of the (i, j) entry, to the lower-tri vector indices.
//      When (i > j), the entry is on the lower triangular part, return the normal mat2vec((i - 1) + (j - 1) * dim).
//      When (i < j), the entry is on the upper triangular part, return the negation of mat2vec((j - 1) + (i - 1) * dim).
//      When (i == j), the entry on diagonal and sometime does not carry info. Return dim * (dim - 1) / 2 + i at the end of vec that may not exist.
// (2) vec2low, size: (dim * (dim - 1) / 2) + dim, it maps the lower-tri vector indices to the column major matrix indices in the lower triangular part.
// (3) vec2upp, size: (dim * (dim - 1) / 2) + dim, it maps the lower-tri vector indices to the column major matrix indices of the upper triangular part.
// (4) In the case when the diagonal carries no info, i.e., strict_lower = true, the diagonals of mat always map to the last part of the vec, i.e.,
//      mat2vec[i + dim * i] = (dim * (dim - 1) / 2) + i
//      vec2low[(dim * (dim - 1) / 2) + i] = i + dim * i
//      vec2upp[(dim * (dim - 1) / 2) + i] = i + dim * i
// Note that (2) and (3) together essentially store the inverse map of (1).

void strict_lower_blk_traversal(INTE_TYPE *, INTE_TYPE *, INTE_TYPE *, INTE_TYPE);
void strict_lower_col_traversal(INTE_TYPE *, INTE_TYPE *, INTE_TYPE *, INTE_TYPE);

class LowerTraversal
{
public:
    INTE_TYPE *m2v;
    INTE_TYPE *v2l;
    INTE_TYPE *v2u;
    INTE_TYPE d;
    INTE_TYPE l;
    BOOL_TYPE sl;

    LowerTraversal() : m2v(nullptr), v2l(nullptr), v2u(nullptr), d(0), l(0), sl(false){};
    LowerTraversal(INTE_TYPE dim, BOOL_TYPE strict_lower, void (*traversal)(INTE_TYPE *, INTE_TYPE *, INTE_TYPE *, INTE_TYPE))
    {
        d = dim;
        sl = strict_lower;
        l = sl ? ((d - 1) * d) / 2 : ((d + 1) * d) / 2;
        m2v = new INTE_TYPE[d * d];
        v2l = new INTE_TYPE[((d + 1) * d) / 2];
        v2u = new INTE_TYPE[((d + 1) * d) / 2];
        traversal(m2v, v2l, v2u, dim);
    };
    LowerTraversal(const LowerTraversal &src) : d(src.d), l(src.l)
    {
        m2v = new INTE_TYPE[d * d];
        v2l = new INTE_TYPE[l];
        v2u = new INTE_TYPE[l];

        memcpy(this->m2v, src.m2v, sizeof(INTE_TYPE) * d * d);
        memcpy(this->v2l, src.v2l, sizeof(INTE_TYPE) * l);
        memcpy(this->v2u, src.v2u, sizeof(INTE_TYPE) * l);
    };
    void swap(LowerTraversal &src)
    {
        using std::swap;
        swap(this->m2v, src.m2v);
        swap(this->v2l, src.v2l);
        swap(this->v2u, src.v2u);
        swap(this->d, src.d);
        swap(this->l, src.l);
    };
    LowerTraversal &operator=(const LowerTraversal &rhs)
    {
        LowerTraversal temp(rhs);
        swap(temp);
        return (*this);
    }
    ~LowerTraversal()
    {
        if (m2v)
            delete[] m2v;
        m2v = nullptr;
        if (v2l)
            delete[] v2l;
        v2l = nullptr;
        if (v2u)
            delete[] v2u;
        v2u = nullptr;
    }
};
