// ================================================
//   mex_sblas_logm.cpp
//   MATLAB usage:
//      X = mex_sblas_logm(Q)
//      [X, R, a] = mex_sblas_logm(Q)
//      X = mex_sblas_logm(Q, e)
//      [X, R, a] = mex_sblas_logm(Q, e)
// ================================================

#include "mex.h"
// #define MATLAB_MEX_BUILD

#include "../src/blasType_mex.hpp"
#include "../src/matOp.hpp"
#include "../src/skewSchFac.hpp"

#include <cstring>
#include <cmath>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // ------------------------------------------------------------
    // Input checking
    // ------------------------------------------------------------
    if (nrhs < 1 || nrhs > 2)
        mexErrMsgIdAndTxt("skewblas:logm:nrhs",
                          "Usage: X = mex_sblas_logm(Q) or X = mex_sblas_logm(Q, e)");

    if (nlhs > 3)
        mexErrMsgIdAndTxt("skewblas:logm:nlhs",
                          "Too many output arguments.");

    const mxArray *Q_in = prhs[0];
    if (!mxIsDouble(Q_in) || mxIsComplex(Q_in))
        mexErrMsgIdAndTxt("skewblas:logm:typeQ",
                          "Q must be a real double matrix.");

    mwSize n = mxGetM(Q_in);
    if (mxGetN(Q_in) != n)
        mexErrMsgIdAndTxt("skewblas:logm:squareQ",
                          "Q must be square.");

    const mwSize m = n / 2;   // floor(n/2)

    bool use_shift = (nrhs == 2);
    const mxArray *E_in = nullptr;

    if (use_shift)
    {
        E_in = prhs[1];

        if (mxIsComplex(E_in))
            mexErrMsgIdAndTxt("skewblas:logm:typeE",
                              "Shift vector e must be real.");

        if (mxGetNumberOfElements(E_in) != m)
            mexErrMsgIdAndTxt("skewblas:logm:sizeE",
                              "Shift vector e must have length floor(n/2).");
    }

    const double *Q = mxGetPr(Q_in);

    // ------------------------------------------------------------
    // Allocate outputs
    // ------------------------------------------------------------
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);     // X
    double *X = mxGetPr(plhs[0]);

    if (nlhs >= 2)
        plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL); // R

    if (nlhs >= 3)
        plhs[2] = mxCreateDoubleMatrix(m, 1, mxREAL); // a

    // ------------------------------------------------------------
    // Modifiable copy of Q
    // ------------------------------------------------------------
    ColMat<double> tmpQ(n, n);
    std::memcpy(tmpQ.v, Q, sizeof(double) * n * n);

    // ------------------------------------------------------------
    // Schur factorization of Q
    // ------------------------------------------------------------
    SkewSchurFactor ssf((int)n);

    ssf.Factor_SpecOrth(tmpQ.v, (int)n);
    ssf.Explict_Vector();

    // ------------------------------------------------------------
    // Case 1: principal logarithm
    // ------------------------------------------------------------
    if (!use_shift)
    {
        ssf.Principal_Logarithm(X, (int)n);
    }
    // ------------------------------------------------------------
    // Case 2: canonical alignment + prescribed 2pi shifts
    // ------------------------------------------------------------
    else
    {
        const double two_pi = 2.0 * M_PI;

        ssf.Canonical_SchurAngular();

        if (mxIsInt32(E_in))
        {
            const int32_T *e_in = static_cast<const int32_T *>(mxGetData(E_in));
            for (mwSize k = 0; k < m; ++k)
            {
                ssf.E[k] = static_cast<INTE_TYPE>(e_in[k]);
                ssf.A[k] += two_pi * static_cast<double>(ssf.E[k]);
            }
        }
        else if (mxIsInt64(E_in))
        {
            const int64_T *e_in = static_cast<const int64_T *>(mxGetData(E_in));
            for (mwSize k = 0; k < m; ++k)
            {
                ssf.E[k] = static_cast<INTE_TYPE>(e_in[k]);
                ssf.A[k] += two_pi * static_cast<double>(ssf.E[k]);
            }
        }
        else if (mxIsDouble(E_in))
        {
            const double *e_in = mxGetPr(E_in);
            for (mwSize k = 0; k < m; ++k)
            {
                double ek = e_in[k];
                double erk = std::round(ek);
                if (std::abs(ek - erk) > 1e-12)
                    mexErrMsgIdAndTxt("skewblas:logm:nonintegerE",
                                      "Shift vector e must contain integers.");

                ssf.E[k] = static_cast<INTE_TYPE>(erk);
                ssf.A[k] += two_pi * static_cast<double>(ssf.E[k]);
            }
        }
        else
        {
            mexErrMsgIdAndTxt("skewblas:logm:typeE",
                              "Shift vector e must be int32, int64, or double-valued integers.");
        }

        ssf.GetSkewSymm(X, (int)n);
    }

    skewl2m(X, (int)n, (int)n);

    // ------------------------------------------------------------
    // Return Schur information if requested
    // ------------------------------------------------------------
    if (nlhs >= 2)
        std::memcpy(mxGetPr(plhs[1]), ssf.R.v, sizeof(double) * n * n);

    if (nlhs >= 3)
        std::memcpy(mxGetPr(plhs[2]), ssf.A.v, sizeof(double) * m);
}