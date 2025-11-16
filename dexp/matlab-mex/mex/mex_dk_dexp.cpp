#include "mex.h"
#include "matrix.h"
#include <string>
#include <cstring>

extern "C"
{
    void DexpDalKreinForward(double *RealMatY, double *RealMatX,
                             double *CmpxMatEigVec, double *CmpxParaForward, int n);

    void DexpDalKreinInverse(double *RealMatY, double *RealMatX,
                             double *CmpxMatEigVec, double *CmpxParaInverse, int n);
}

static void fail(const char *msg)
{
    mexErrMsgIdAndTxt("skewblas:dkdexp", msg);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /*
        Usage:
            Y = mex_dk_dexp(EigVec, PF, PI, X, dir)

        where:
            - EigVec : 2*n*n real vector (interleaved complex eigenvectors)
            - PF, PI : 2*n*n real vectors (interleaved DK parameters)
            - X      : n×n real matrix
            - dir    : 'fwd' or 'inv'
    */

    if (nrhs != 5)
        fail("Usage: Y = mex_dk_dexp(EigVec, PF, PI, X, dir)");

    /* ============================================================
       1) Read interleaved EigVec (2*n*n real vector)
       ============================================================ */
    const mxArray *EigVec_mx = prhs[0];
    if (!mxIsDouble(EigVec_mx) || mxIsComplex(EigVec_mx))
        fail("EigVec must be a real 2*n*n vector (interleaved).");

    mwSize lenV = mxGetNumberOfElements(EigVec_mx);
    int n = static_cast<int>(sqrt(lenV / 2.0) + 1e-12);

    if (2 * n * n != lenV)
        fail("EigVec length must be 2*n*n (real interleaved complex data).");

    double *Vint = mxGetPr(EigVec_mx);

    /* ============================================================
       2) Read PF, PI (real interleaved 2*n*n vectors)
       ============================================================ */
    const mxArray *PF_mx = prhs[1];
    const mxArray *PI_mx = prhs[2];

    if (mxGetNumberOfElements(PF_mx) != 2 * n * n ||
        mxGetNumberOfElements(PI_mx) != 2 * n * n)
        fail("PF and PI must be 2*n*n real vectors (interleaved).");

    double *PF = mxGetPr(PF_mx);
    double *PI = mxGetPr(PI_mx);

    /* ============================================================
       3) Read X (real n×n)
       ============================================================ */
    const mxArray *X_mx = prhs[3];

    if (!mxIsDouble(X_mx) || mxIsComplex(X_mx))
        fail("X must be real double.");

    if ((int)mxGetM(X_mx) != n || (int)mxGetN(X_mx) != n)
        fail("X must be n×n.");

    double *X = mxGetPr(X_mx);

    /* ============================================================
       4) Read direction flag ('fwd' or 'inv')
       ============================================================ */
    char dirbuf[16];
    mxGetString(prhs[4], dirbuf, sizeof(dirbuf));
    std::string dir(dirbuf);

    bool is_forward;
    if (dir == "fwd")
        is_forward = true;
    else if (dir == "inv")
        is_forward = false;
    else
        fail("dir must be 'fwd' or 'inv'.");

    /* ============================================================
       5) Allocate output Y (real n×n)
       ============================================================ */
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    double *Y = mxGetPr(plhs[0]);

    /* ============================================================
       6) Dispatch backend call
       ============================================================ */
    if (is_forward)
    {
        DexpDalKreinForward(Y, X, Vint, PF, n);
    }
    else
    {
        DexpDalKreinInverse(Y, X, Vint, PI, n);
    }
}
