// ================================================
//   mex_sblas_expm.cpp
//   MATLAB usage:
//      Q = mex_sblas_expm(S)
// ================================================

#include "mex.h"
// #define MATLAB_MEX_BUILD

#include "../src/blasType_mex.hpp"
#include "../src/matOp.hpp"
#include "../src/skewSchFac.hpp"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // ------------------------------------------------------------
    // Check inputs
    // ------------------------------------------------------------
    if (nrhs != 1)
        mexErrMsgIdAndTxt("skewblas:expm:nrhs",
                          "Usage: Q = mex_sblas_expm(S)");

    const mxArray *S_in = prhs[0];
    if (!mxIsDouble(S_in) || mxIsComplex(S_in))
        mexErrMsgIdAndTxt("skewblas:expm:type",
                          "Input S must be a real double matrix.");

    mwSize n = mxGetM(S_in);
    if (mxGetN(S_in) != n)
        mexErrMsgIdAndTxt("skewblas:expm:square",
                          "S must be square.");

    const double *S = mxGetPr(S_in);

    // ------------------------------------------------------------
    // Allocate output Q
    // ------------------------------------------------------------
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    double *Q = mxGetPr(plhs[0]);

    // ------------------------------------------------------------
    // Perform expm(S) via C++ class SkewSchurFactor
    // ------------------------------------------------------------
    SkewSchurFactor ssf((int)n);

    // Copy S into the internal buffer
    memcpy(ssf.H.MatH.v, S, n * n * sizeof(double));

    ssf.SchurAngular_SkewSymm();
    ssf.Explict_Vector();
    ssf.Exponential(Q, (int)n);
}
