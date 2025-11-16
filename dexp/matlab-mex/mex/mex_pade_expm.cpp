// ================================================
//   mex_pade_expm.cpp
//   MATLAB usage:
//      Q = mex_pade_expm(A, p)
// ================================================

#include "mex.h"
// #define MATLAB_MEX_BUILD

#include "../src/blasType_mex.hpp"
#include "../src/matOp.hpp"
#include "../src/expmPade.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // ------------------------------------------------------------
    // Input checking
    // ------------------------------------------------------------
    if (nrhs < 1 || nrhs > 2)
        mexErrMsgIdAndTxt("skewblas:pade:nrhs",
                          "Usage: Q = mex_pade_expm(A, p)");

    const mxArray *A_in = prhs[0];
    if (!mxIsDouble(A_in) || mxIsComplex(A_in))
        mexErrMsgIdAndTxt("skewblas:pade:type",
                          "A must be a real double matrix.");

    mwSize n = mxGetM(A_in);
    if (mxGetN(A_in) != n)
        mexErrMsgIdAndTxt("skewblas:pade:square",
                          "A must be square.");

    int p = 7; // default Pade order
    if (nrhs == 2)
        p = (int)mxGetScalar(prhs[1]);

    const double *A = mxGetPr(A_in);

    // ------------------------------------------------------------
    // Allocate output Q
    // ------------------------------------------------------------
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    double *Q = mxGetPr(plhs[0]);

    // ------------------------------------------------------------
    // Compute exp(A) using Padé approximation
    // ------------------------------------------------------------
    expmPadeApprox pade((int)n);

    // internally pade.Expm(Q, n, A, n, p)
    // expects a modifiable copy of A
    ColMat<double> tmpA(n, n);
    memcpy(tmpA.v, A, n * n * sizeof(double));

    pade.Expm(Q, (int)n, tmpA.v, (int)n, p);
}
