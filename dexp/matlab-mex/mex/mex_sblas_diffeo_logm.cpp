// ================================================
//   mex_sblas_diffeo_logm.cpp
//   MATLAB usage:
//      [X, R, a, e, V, b] = mex_sblas_diffeo_logm(Q, S)
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
        mexErrMsgIdAndTxt("skewblas:diffeo_logm:nrhs",
                          "Usage: [X, R, a, e, V, b] = mex_sblas_diffeo_logm(Q, S)");

    if (nlhs > 6)
        mexErrMsgIdAndTxt("skewblas:diffeo_logm:nlhs",
                          "Too many output arguments.");

    const mxArray *Q_in = prhs[0];
    const mxArray *S_in = prhs[1];

    if (!mxIsDouble(Q_in) || mxIsComplex(Q_in))
        mexErrMsgIdAndTxt("skewblas:diffeo_logm:typeQ",
                          "Q must be a real double matrix.");

    if (!mxIsDouble(S_in) || mxIsComplex(S_in))
        mexErrMsgIdAndTxt("skewblas:diffeo_log:typeS",
                          "S must be a real double matrix.");

    mwSize n = mxGetM(Q_in);
    if (mxGetN(Q_in) != n)
        mexErrMsgIdAndTxt("skewblas:diffeo_log:squareQ",
                          "Q must be square.");

    if (mxGetM(S_in) != n || mxGetN(S_in) != n)
        mexErrMsgIdAndTxt("skewblas:diffeo_log:size",
                          "Q and S must have the same size.");

    const double *Q = mxGetPr(Q_in);
    const double *S = mxGetPr(S_in);

    // ------------------------------------------------------------
    // Allocate outputs
    // ------------------------------------------------------------
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);          // X
    double *X = mxGetPr(plhs[0]);

    if (nlhs >= 2)
        plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);      // R

    if (nlhs >= 3)
        plhs[2] = mxCreateDoubleMatrix(n / 2, 1, mxREAL);  // a

    if (nlhs >= 4)
        plhs[3] = mxCreateNumericMatrix(n / 2, 1, mxINT32_CLASS, mxREAL); // e

    if (nlhs >= 5)
        plhs[4] = mxCreateDoubleMatrix(n, n, mxREAL);      // V

    if (nlhs >= 6)
        plhs[5] = mxCreateDoubleMatrix(n / 2, 1, mxREAL);  // b

    // ------------------------------------------------------------
    // Modifiable copies
    // ------------------------------------------------------------
    ColMat<double> tmpQ(n, n);
    ColMat<double> tmpS(n, n);

    std::memcpy(tmpQ.v, Q, sizeof(double) * n * n);
    std::memcpy(tmpS.v, S, sizeof(double) * n * n);

    // ------------------------------------------------------------
    // Build Schur factors
    // ------------------------------------------------------------
    SkewSchurFactor ssfQ((int)n);
    SkewSchurFactor ssfS((int)n);

    ssfQ.Factor_SpecOrth(tmpQ.v, (int)n);
    ssfQ.Explict_Vector();

    ssfS.Factor_SkewSymm(tmpS.v, (int)n);
    ssfS.Explict_Vector();

    // ------------------------------------------------------------
    // Canonical alignment required before calling Diffeomorphic_Logarithm
    // ------------------------------------------------------------
    ssfQ.Canonical_SchurAngular();
    ssfS.Canonical_SchurAngular();

    // ------------------------------------------------------------
    // Compute diffeomorphic logarithm
    // ------------------------------------------------------------
    ssfQ.Diffeomorphic_Logarithm(X, (int)n, tmpS.v, (int)n, ssfS);

    skewl2m(X, (int)n, (int)n);

    // ------------------------------------------------------------
    // Return outputs
    // ------------------------------------------------------------
    if (nlhs >= 2)
        std::memcpy(mxGetPr(plhs[1]), ssfQ.R.v, sizeof(double) * n * n);

    if (nlhs >= 3)
        std::memcpy(mxGetPr(plhs[2]), ssfQ.A.v, sizeof(double) * (n / 2));

    if (nlhs >= 4)
    {
        int32_T *e_out = static_cast<int32_T *>(mxGetData(plhs[3]));
        for (mwSize i = 0; i < n / 2; ++i)
            e_out[i] = static_cast<int32_T>(ssfQ.E[i]);
    }

    if (nlhs >= 5)
        std::memcpy(mxGetPr(plhs[4]), ssfS.R.v, sizeof(double) * n * n);

    if (nlhs >= 6)
        std::memcpy(mxGetPr(plhs[5]), ssfS.A.v, sizeof(double) * (n / 2));
}