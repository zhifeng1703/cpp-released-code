// ================================================
//   mex_pade_expm.cpp
//   MATLAB usage:
//      Q = mex_pade_expm(A, p)
//      [Q, info] = mex_pade_expm(A, p)
// ================================================

#include "mex.h"
// #define MATLAB_MEX_BUILD

#include "../src/blasType_mex.hpp"
#include "../src/matOp.hpp"
#include "../src/expmPade.hpp"

static mxArray *export_dyadic_cache(const expmPadeApprox &pade, mwSize n)
{
    // Export in MATLAB-friendly order:
    //   dyadic{1}   = exp(A)
    //   dyadic{2}   = exp(A/2)
    //   ...
    //   dyadic{s+1} = exp(A/2^s)
    mwSize L = static_cast<mwSize>(pade.s + 1);
    mxArray *cell = mxCreateCellMatrix(L, 1);

    for (mwIndex k = 0; k < L; ++k)
    {
        // Internal storage:
        //   VecMatQ[0] = exp(A/2^s)
        //   VecMatQ[s] = exp(A)
        mwIndex internal_idx = static_cast<mwIndex>(pade.s) - k;

        mxArray *Mk = mxCreateDoubleMatrix(n, n, mxREAL);
        double *dst = mxGetPr(Mk);
        memcpy(dst, pade.VecMatQ[internal_idx]->v, n * n * sizeof(double));
        mxSetCell(cell, k, Mk);
    }

    return cell;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // ------------------------------------------------------------
    // Input checking
    // ------------------------------------------------------------
    if (nrhs < 1 || nrhs > 2)
        mexErrMsgIdAndTxt("skewblas:pade:nrhs",
                          "Usage: Q = mex_pade_expm(A, p) or [Q, info] = mex_pade_expm(A, p)");

    if (nlhs > 2)
        mexErrMsgIdAndTxt("skewblas:pade:nlhs",
                          "Too many output arguments.");

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
        p = static_cast<int>(mxGetScalar(prhs[1]));

    const double *A = mxGetPr(A_in);

    // ------------------------------------------------------------
    // Allocate output Q
    // ------------------------------------------------------------
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    double *Q = mxGetPr(plhs[0]);

    // ------------------------------------------------------------
    // Compute exp(A) using Padé approximation
    // ------------------------------------------------------------
    expmPadeApprox pade(static_cast<INTE_TYPE>(n));

    // pade.Expm expects a modifiable copy of A
    ColMat<double> tmpA(static_cast<INTE_TYPE>(n), static_cast<INTE_TYPE>(n));
    memcpy(tmpA.v, A, n * n * sizeof(double));

    pade.Expm(Q, static_cast<INTE_TYPE>(n), tmpA.v, static_cast<INTE_TYPE>(n), static_cast<INTE_TYPE>(p));

    // ------------------------------------------------------------
    // Optional second output: expose internal cache and metadata
    // ------------------------------------------------------------
    if (nlhs >= 2)
    {
        const char *fields[] = {"dyadic", "depth", "order", "selected_order"};
        mxArray *info = mxCreateStructMatrix(1, 1, 4, fields);

        // dyadic{1} = exp(A), dyadic{2} = exp(A/2), ...
        mxSetField(info, 0, "dyadic", export_dyadic_cache(pade, n));

        // Actual scaling depth used
        mxSetField(info, 0, "depth", mxCreateDoubleScalar(static_cast<double>(pade.s)));

        // Requested order from input
        mxSetField(info, 0, "order", mxCreateDoubleScalar(static_cast<double>(p)));

        // Actual order selected internally (useful if you later call the automatic-order path)
        mxSetField(info, 0, "selected_order", mxCreateDoubleScalar(static_cast<double>(pade.m)));

        plhs[1] = info;
    }
}