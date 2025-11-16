// ============================================================================
//   blasType_mex.hpp  (Corrected for MATLAB R2024a, MSVC 2022)
//   BLAS/LAPACK shim to replace MKL when MATLAB_MEX_BUILD is defined
//
//   - Uses MWBLAS & MWLAPACK (Fortran symbols in MATLAB)
//   - Corrects signatures for zgemm, dbdsdc, dsytrd, dorgtr, dormtr
//   - Compatible with your skewblas code without modification
//
// ============================================================================

#ifndef BLASTYPE_MEX_HPP
#define BLASTYPE_MEX_HPP

#include <cmath>
#include <cstddef>
#include <vector> // IMPORTANT — needed for workspaces
#include "mex.h"

extern "C"
{
#include "blas.h"
#include "lapack.h"
}

// ---------------------------------------------------------------------
// Basic types
// ---------------------------------------------------------------------
#define BOOL_TYPE bool
#define REAL_TYPE double
#define INTE_TYPE ptrdiff_t // MATLAB BLAS/LAPACK want ptrdiff_t
#define CHAR_TYPE char

typedef struct
{
    double real;
    double imag;
} CMPX_TYPE;

// ---------------------------------------------------------------------
// Complex arithmetic
// ---------------------------------------------------------------------

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

inline void expm1divCMPX(CMPX_TYPE &res, const CMPX_TYPE &z)
{
    const double tol = 1e-8;
    double r2 = z.real * z.real + z.imag * z.imag;
    if (r2 < tol * tol)
    {
        // series expansion up to z^2
        CMPX_TYPE term;
        assignCMPX(res, 1.0 + 0.5 * z.real, 0.5 * z.imag);
        term.real = (z.real * z.real - z.imag * z.imag) / 6.0;
        term.imag = (2.0 * z.real * z.imag) / 6.0;
        additiCMPX(res, res, term);
    }
    else
    {
        CMPX_TYPE ez;
        exponeCMPX(ez, z);
        // ez - 1
        ez.real -= 1.0;
        divideCMPX(res, ez, z);
    }
}

// ---------------------------------------------------------------------
// CBLAS enums
// ---------------------------------------------------------------------
enum CBLAS_LAYOUT
{
    CblasRowMajor = 101,
    CblasColMajor = 102
};
enum CBLAS_TRANSPOSE
{
    CblasNoTrans = 111,
    CblasTrans = 112,
    CblasConjTrans = 113
};

#ifndef LAPACK_COL_MAJOR
#define LAPACK_COL_MAJOR 102
#endif

inline CHAR_TYPE trans_to_char(CBLAS_TRANSPOSE t)
{
    switch (t)
    {
    case CblasNoTrans:
        return 'N';
    case CblasTrans:
        return 'T';
    case CblasConjTrans:
        return 'C';
    }
    return 'N';
}

// ============================================================================
// LEVEL-1 BLAS wrappers
// ============================================================================
inline double cblas_dasum(INTE_TYPE n, const double *x, INTE_TYPE incx)
{
    return dasum(&n, x, &incx);
}

inline void cblas_daxpy(INTE_TYPE n, double alpha,
                        const double *x, INTE_TYPE incx,
                        double *y, INTE_TYPE incy)
{
    daxpy(&n, &alpha, x, &incx, y, &incy);
}

// MKL extension: implement manually
inline void cblas_daxpby(INTE_TYPE n, double alpha,
                         const double *x, INTE_TYPE incx,
                         double beta,
                         double *y, INTE_TYPE incy)
{
    INTE_TYPE ix = 0, iy = 0;
    for (INTE_TYPE i = 0; i < n; ++i)
    {
        y[iy] = alpha * x[ix] + beta * y[iy];
        ix += incx;
        iy += incy;
    }
}

inline void cblas_dscal(INTE_TYPE n, double alpha, double *x, INTE_TYPE incx)
{
    dscal(&n, &alpha, x, &incx);
}

inline void cblas_dcopy(INTE_TYPE n, const double *x, INTE_TYPE incx,
                        double *y, INTE_TYPE incy)
{
    dcopy(&n, x, &incx, y, &incy);
}

inline double cblas_dnrm2(INTE_TYPE n, const double *x, INTE_TYPE incx)
{
    return dnrm2(&n, x, &incx);
}

inline double cblas_ddot(INTE_TYPE n, const double *x, INTE_TYPE incx,
                         const double *y, INTE_TYPE incy)
{
    return ddot(&n, x, &incx, y, &incy);
}

// ============================================================================
// LEVEL-2 BLAS wrappers
// ============================================================================
inline void cblas_dgemv(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE trans,
                        INTE_TYPE m, INTE_TYPE n, double alpha,
                        const double *A, INTE_TYPE lda,
                        const double *x, INTE_TYPE incx,
                        double beta,
                        double *y, INTE_TYPE incy)
{
    if (layout != CblasColMajor)
        mexErrMsgIdAndTxt("skewblas:dgemv", "Only column-major supported");

    CHAR_TYPE ta = trans_to_char(trans);
    dgemv(&ta, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}

inline void cblas_dger(CBLAS_LAYOUT layout,
                       INTE_TYPE m, INTE_TYPE n,
                       double alpha,
                       const double *x, INTE_TYPE incx,
                       const double *y, INTE_TYPE incy,
                       double *A, INTE_TYPE lda)
{
    if (layout != CblasColMajor)
        mexErrMsgIdAndTxt("skewblas:dger", "Only column-major supported");

    dger(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
}

// ============================================================================
// LEVEL-3 BLAS wrappers
// ============================================================================

inline void cblas_dgemm(CBLAS_LAYOUT layout,
                        CBLAS_TRANSPOSE transA,
                        CBLAS_TRANSPOSE transB,
                        INTE_TYPE m, INTE_TYPE n, INTE_TYPE k,
                        double alpha,
                        const double *A, INTE_TYPE lda,
                        const double *B, INTE_TYPE ldb,
                        double beta,
                        double *C, INTE_TYPE ldc)
{
    if (layout != CblasColMajor)
        mexErrMsgIdAndTxt("skewblas:dgemm", "Only column-major supported");

    CHAR_TYPE ta = trans_to_char(transA);
    CHAR_TYPE tb = trans_to_char(transB);

    dgemm(&ta, &tb, &m, &n, &k,
          &alpha,
          A, &lda,
          B, &ldb,
          &beta,
          C, &ldc);
}

// ---- COMPLEX GEMM fix (MATLAB expects double* interleaved) ----
inline void cblas_zgemm(CBLAS_LAYOUT layout,
                        CBLAS_TRANSPOSE transA,
                        CBLAS_TRANSPOSE transB,
                        INTE_TYPE m, INTE_TYPE n, INTE_TYPE k,
                        const CMPX_TYPE *alpha,
                        const CMPX_TYPE *A, INTE_TYPE lda,
                        const CMPX_TYPE *B, INTE_TYPE ldb,
                        const CMPX_TYPE *beta,
                        CMPX_TYPE *C, INTE_TYPE ldc)
{
    if (layout != CblasColMajor)
        mexErrMsgIdAndTxt("skewblas:zgemm", "Only column-major supported");

    CHAR_TYPE ta = trans_to_char(transA);
    CHAR_TYPE tb = trans_to_char(transB);

    const double *alpha_raw = reinterpret_cast<const double *>(alpha);
    const double *beta_raw = reinterpret_cast<const double *>(beta);
    const double *A_raw = reinterpret_cast<const double *>(A);
    const double *B_raw = reinterpret_cast<const double *>(B);
    double *C_raw = reinterpret_cast<double *>(C);

    zgemm(&ta, &tb, &m, &n, &k,
          alpha_raw,
          A_raw, &lda,
          B_raw, &ldb,
          beta_raw,
          C_raw, &ldc);
}

// ============================================================================
// LAPACK wrappers
// ============================================================================

// ---- dgetrf ----
inline INTE_TYPE LAPACKE_dgetrf(int layout,
                                INTE_TYPE m, INTE_TYPE n,
                                double *A, INTE_TYPE lda,
                                INTE_TYPE *ipiv)
{
    if (layout != LAPACK_COL_MAJOR)
        mexErrMsgIdAndTxt("skewblas:dgetrf", "Only column-major supported");

    ptrdiff_t M = m, N = n, LDA = lda;
    ptrdiff_t INFO = 0;

    dgetrf(&M, &N, A, &LDA, ipiv, &INFO);

    return INFO;
}

// ---- dgetrs ----
inline INTE_TYPE LAPACKE_dgetrs(int layout, char trans,
                                INTE_TYPE n, INTE_TYPE nrhs,
                                double *A, INTE_TYPE lda,
                                INTE_TYPE *ipiv,
                                double *B, INTE_TYPE ldb)
{
    if (layout != LAPACK_COL_MAJOR)
        mexErrMsgIdAndTxt("skewblas:dgetrs", "Only column-major supported");

    ptrdiff_t N = n, NRHS = nrhs, LDA = lda, LDB = ldb;
    ptrdiff_t INFO = 0;

    dgetrs(&trans, &N, &NRHS, A, &LDA, ipiv, B, &LDB, &INFO);
    return INFO;
}

// ---- dbdsdc (correct signature) ----
inline INTE_TYPE LAPACKE_dbdsdc(int layout, char uplo, char compq,
                                INTE_TYPE n,
                                double *D, double *E,
                                double *U, INTE_TYPE ldu,
                                double *VT, INTE_TYPE ldvt,
                                double *Q, INTE_TYPE *IQ)
{
    if (layout != LAPACK_COL_MAJOR)
        mexErrMsgIdAndTxt("skewblas:dbdsdc", "Only column-major supported");

    ptrdiff_t N = n, LDU = ldu, LDVT = ldvt;
    ptrdiff_t INFO = 0;

    ptrdiff_t LWORK = 0;
    if (compq == 'N')
        LWORK = 4 * N;
    else if (compq == 'P')
        LWORK = 6 * N;
    else
        LWORK = 3 * N * N + 4 * N;
    ptrdiff_t LIWORK = 8 * N;

    // std::vector<double> work(LWORK);
    // std::vector<ptrdiff_t> iwork(LIWORK);

    double *work = new double[LWORK];
    ptrdiff_t *iwork = new ptrdiff_t[LIWORK];

    dbdsdc(&uplo, &compq, &N,
           D, E,
           U, &LDU,
           VT, &LDVT,
           Q, IQ,
           work, iwork,
           &INFO);

    delete[] work, iwork;

    return INFO;
}

// ---- dsytrd ----
inline INTE_TYPE LAPACKE_dsytrd(int layout, char uplo,
                                INTE_TYPE n,
                                double *A, INTE_TYPE lda,
                                double *d, double *e, double *tau)
{
    if (layout != LAPACK_COL_MAJOR)
        mexErrMsgIdAndTxt("skewblas:dsytrd", "Only column-major supported");

    ptrdiff_t N = n, LDA = lda, INFO = 0, LWORK = -1;
    double workq;

    // Workspace query
    dsytrd(&uplo, &N, A, &LDA, d, e, tau, &workq, &LWORK, &INFO);

    LWORK = (ptrdiff_t)workq;
    std::vector<double> work(LWORK);

    dsytrd(&uplo, &N, A, &LDA, d, e, tau, work.data(), &LWORK, &INFO);

    return INFO;
}

// ---- dorgtr ----
inline INTE_TYPE LAPACKE_dorgtr(int layout, char uplo,
                                INTE_TYPE n,
                                double *A, INTE_TYPE lda,
                                double *tau)
{
    if (layout != LAPACK_COL_MAJOR)
        mexErrMsgIdAndTxt("skewblas:dorgtr", "Only column-major supported");

    ptrdiff_t N = n, LDA = lda, INFO = 0, LWORK = -1;
    double workq;

    dorgtr(&uplo, &N, A, &LDA, tau, &workq, &LWORK, &INFO);

    LWORK = (ptrdiff_t)workq;
    std::vector<double> work(LWORK);

    dorgtr(&uplo, &N, A, &LDA, tau, work.data(), &LWORK, &INFO);

    return INFO;
}

// ---- dormtr ----
inline INTE_TYPE LAPACKE_dormtr(int layout, char side, char uplo,
                                char trans,
                                INTE_TYPE m, INTE_TYPE n,
                                double *A, INTE_TYPE lda,
                                double *tau,
                                double *C, INTE_TYPE ldc)
{
    if (layout != LAPACK_COL_MAJOR)
        mexErrMsgIdAndTxt("skewblas:dormtr", "Only column-major supported");

    ptrdiff_t M = m, N = n, LDA = lda, LDC = ldc;
    ptrdiff_t INFO = 0, LWORK = -1;
    double workq;

    dormtr(&side, &uplo, &trans,
           &M, &N, A, &LDA,
           tau, C, &LDC,
           &workq, &LWORK, &INFO);

    LWORK = (ptrdiff_t)workq;
    std::vector<double> work(LWORK);

    dormtr(&side, &uplo, &trans,
           &M, &N, A, &LDA,
           tau, C, &LDC,
           work.data(), &LWORK, &INFO);

    return INFO;
}

#endif // BLASTYPE_MEX_HPP
