#include "mex.h"
#include "matrix.h"

extern "C"
{
    void DexpPadeForward(double *MatY, double *MatX,
                         double *MatS, int n, int p);
}

static void fail(const char *msg)
{
    mexErrMsgIdAndTxt("skewblas:pade_fwd", msg);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3)
        fail("Usage: Y = mex_pade_dexpfwd(S, X, p)");

    /* ----------------------------------
       Read S
       ---------------------------------- */
    const mxArray *S_mx = prhs[0];
    if (!mxIsDouble(S_mx) || mxIsComplex(S_mx))
        fail("S must be real double.");

    int n = (int)mxGetM(S_mx);
    if (mxGetN(S_mx) != n)
        fail("S must be n×n.");

    double *S = mxGetPr(S_mx);

    /* ----------------------------------
       Read X
       ---------------------------------- */
    const mxArray *X_mx = prhs[1];
    if (!mxIsDouble(X_mx) || mxIsComplex(X_mx))
        fail("X must be real double.");

    if ((int)mxGetM(X_mx) != n || (int)mxGetN(X_mx) != n)
        fail("X must be n×n.");

    double *X = mxGetPr(X_mx);

    /* ----------------------------------
       Read p
       ---------------------------------- */
    const mxArray *p_mx = prhs[2];
    if (!mxIsDouble(p_mx) || mxIsComplex(p_mx) || mxGetNumberOfElements(p_mx) != 1)
        fail("p must be a real scalar.");

    int p = (int)*mxGetPr(p_mx);

    /* ----------------------------------
       Output Y
       ---------------------------------- */
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    double *Y = mxGetPr(plhs[0]);

    /* ----------------------------------
       Call backend
       ---------------------------------- */
    DexpPadeForward(Y, X, S, n, p);
}
