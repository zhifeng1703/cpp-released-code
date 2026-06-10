#include "mex.h"
#include "matrix.h"

extern "C"
{
    void DexpPadeInverse(double *MatY, double *MatX,
                         double *MatR, double *VecA,
                         int dim, int pade_order, int scale_order);
}

static void fail(const char *msg)
{
    mexErrMsgIdAndTxt("skewblas:pade_inv", msg);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 5)
        fail("Usage: Y = mex_pade_dexpinv(R, A, X, p, s)");

    /* ----------------------------------
       Read R
       ---------------------------------- */
    const mxArray *R_mx = prhs[0];
    if (!mxIsDouble(R_mx) || mxIsComplex(R_mx))
        fail("R must be real double.");

    int n = (int)mxGetM(R_mx);
    if (mxGetN(R_mx) != n)
        fail("R must be n×n.");

    double *R = mxGetPr(R_mx);

    /* ----------------------------------
       Read A
       ---------------------------------- */
    const mxArray *A_mx = prhs[1];
    if (!mxIsDouble(A_mx) || mxIsComplex(A_mx))
        fail("A must be real double.");

    if ((int)mxGetNumberOfElements(A_mx) != n / 2)
        fail("A must have n/2 elements.");

    double *A = mxGetPr(A_mx);

    /* ----------------------------------
       Read X
       ---------------------------------- */
    const mxArray *X_mx = prhs[2];
    if (!mxIsDouble(X_mx) || mxIsComplex(X_mx))
        fail("X must be real double.");

    if ((int)mxGetM(X_mx) != n || (int)mxGetN(X_mx) != n)
        fail("X must be n×n.");

    double *X = mxGetPr(X_mx);

    /* ----------------------------------
       p (scalar)
       ---------------------------------- */
    const mxArray *p_mx = prhs[3];
    if (!mxIsDouble(p_mx) || mxIsComplex(p_mx) ||
        mxGetNumberOfElements(p_mx) != 1)
        fail("p must be scalar integer.");

    int p = (int)*mxGetPr(p_mx);

    /* ----------------------------------
    normS (scalar)
    ---------------------------------- */
    const mxArray *s_mx = prhs[4];
    if (!mxIsDouble(s_mx) || mxIsComplex(s_mx) ||
        mxGetNumberOfElements(s_mx) != 1)
        fail("s must be order integer.");

    int s = (int)*mxGetPr(s_mx);

    /* ----------------------------------
       Output Y
       ---------------------------------- */
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    double *Y = mxGetPr(plhs[0]);

    /* ----------------------------------
       Call backend
       ---------------------------------- */

    // mexPrintf("dim %3d, pade order %3d, scale order %3d.\n", n, p, s);

    DexpPadeInverse(Y, X, R, A, n, p, s);
}
