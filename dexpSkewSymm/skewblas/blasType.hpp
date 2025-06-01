#pragma once

#include <cmath>
#include "mkl_types.h"
#include "mkl.h"
// #include "mkl.fi"

#define BOOL_TYPE int
#define REAL_TYPE double
#define INTE_TYPE MKL_INT
#define CMPX_TYPE MKL_Complex16
#define CHAR_TYPE char

inline REAL_TYPE normofCMPX(const CMPX_TYPE &c)
{
    return sqrt(c.real * c.real + c.imag * c.imag);
}
inline void assignCMPX(CMPX_TYPE &c, const REAL_TYPE r, const REAL_TYPE i)
{
    c.real = r;
    c.imag = i;
}
inline void assignCMPX(CMPX_TYPE &lhs, const CMPX_TYPE &rhs)
{
    lhs.real = rhs.real;
    lhs.imag = rhs.imag;
}
inline void additiCMPX(CMPX_TYPE &c, const CMPX_TYPE &a, const CMPX_TYPE &b)
{
    c.real = a.real + b.real;
    c.imag = a.imag + b.imag;
}
inline void substrCMPX(CMPX_TYPE &c, const CMPX_TYPE &a, const CMPX_TYPE &b)
{
    c.real = a.real - b.real;
    c.imag = a.imag - b.imag;
}
inline void multifCMPX(CMPX_TYPE &c, const CMPX_TYPE &a, const CMPX_TYPE &b)
{
    REAL_TYPE r, i;
    r = a.real * b.real - a.imag * b.imag;
    i = a.real * b.imag + b.real * a.imag;
    c.real = r;
    c.imag = i;
}
inline void multifCMPX(CMPX_TYPE &c, REAL_TYPE a, const CMPX_TYPE &b)
{
    c.real = a * b.real;
    c.imag = a * b.imag;
}
inline void divideCMPX(CMPX_TYPE &c, const CMPX_TYPE &a, const CMPX_TYPE &b)
{
    REAL_TYPE bnorm = b.real * b.real + b.imag * b.imag;
    REAL_TYPE r, i;
    r = (a.real * b.real + a.imag * b.imag) / bnorm;
    i = (b.real * a.imag - b.imag * a.real) / bnorm;
    c.real = r;
    c.imag = i;
}
inline void inversCMPX(CMPX_TYPE &y, const CMPX_TYPE &x)
{
    REAL_TYPE xnorm = x.real * x.real + x.imag * x.imag;
    REAL_TYPE r, i;

    r = x.real / xnorm;
    i = -x.imag / xnorm;
    y.real = r;
    y.imag = i;
}
inline void exponeCMPX(CMPX_TYPE &y, const CMPX_TYPE &x)
{
    REAL_TYPE er = exp(x.real);
    y.real = er * cos(x.imag);
    y.imag = er * sin(x.imag);
}
// std::ostream &operator<<(std::ostream &out, CMPX_TYPE c)
// {
//     out << c.real << " + i * " << c.imag;
//     return out;
// }