// ================================================
//   mex_sblas_spdiff_logm.cpp
//   MATLAB usage:
//      [X, R, a] = mex_sblas_spdiff_logm(Q, S)
// ================================================

#include "mex.h"
// #define MATLAB_MEX_BUILD

#include "../src/blasType_mex.hpp"
#include "../src/matOp.hpp"
#include "../src/skewSchFac.hpp"
#include <cstring>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // ------------------------------------------------------------
    // Input checking
    // ------------------------------------------------------------
    if (nrhs != 2)
        mexErrMsgIdAndTxt("skewblas:spdiff_logm:nrhs",
                          "Usage: [X, R, a] = mex_sblas_specialdiffeo_log(Q, S)");

    if (nlhs > 3)
        mexErrMsgIdAndTxt("skewblas:spdiff_logm:nlhs",
                          "Too many output arguments.");

    const mxArray *Q_in = prhs[0];
    const mxArray *S_in = prhs[1];

    if (!mxIsDouble(Q_in) || mxIsComplex(Q_in))
        mexErrMsgIdAndTxt("skewblas:spdiff_logm:typeQ",
                          "Q must be a real double matrix.");

    if (!mxIsDouble(S_in) || mxIsComplex(S_in))
        mexErrMsgIdAndTxt("skewblas:spdiff_logm:typeS",
                          "S must be a real double matrix.");

    mwSize n = mxGetM(Q_in);
    if (mxGetN(Q_in) != n)
        mexErrMsgIdAndTxt("skewblas:spdiff_logm:squareQ",
                          "Q must be square.");

    if (mxGetM(S_in) != n || mxGetN(S_in) != n)
        mexErrMsgIdAndTxt("skewblas:spdiff_logm:size",
                          "Q and S must have the same size.");

    const double *Q = mxGetPr(Q_in);
    const double *S = mxGetPr(S_in);

    // ------------------------------------------------------------
    // Allocate outputs
    // ------------------------------------------------------------
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);     // X
    double *X = mxGetPr(plhs[0]);

    if (nlhs >= 2)
        plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL); // R

    if (nlhs >= 3)
        plhs[2] = mxCreateDoubleMatrix(n / 2, 1, mxREAL); // a

    // ------------------------------------------------------------
    // Modifiable copies
    // ------------------------------------------------------------
    ColMat<double> tmpQ(n, n);
    ColMat<double> tmpS(n, n);

    std::memcpy(tmpQ.v, Q, sizeof(double) * n * n);
    std::memcpy(tmpS.v, S, sizeof(double) * n * n);

    // ------------------------------------------------------------
    // Factor Q and compute special diffeomorphic logarithm
    // ------------------------------------------------------------
    SkewSchurFactor ssfQ((int)n);

    ssfQ.Factor_SpecOrth(tmpQ.v, (int)n);
    ssfQ.Explict_Vector();
    ssfQ.Canonical_SchurAngular();

    ssfQ.SpecialDiffeo_Logarithm(X, (int)n, tmpS.v, (int)n);
    skewl2m(X, (int)n, (int)n);


    // ------------------------------------------------------------
    // Return R and a
    // ------------------------------------------------------------
    if (nlhs >= 2)
        std::memcpy(mxGetPr(plhs[1]), ssfQ.R.v, sizeof(double) * n * n);

    if (nlhs >= 3)
        std::memcpy(mxGetPr(plhs[2]), ssfQ.A.v, sizeof(double) * (n / 2));
}