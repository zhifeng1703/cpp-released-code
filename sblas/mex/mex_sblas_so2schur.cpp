#include "mex.h"
#include "matrix.h"

// C interface provided by skewblas_mex.cpp (compiled from ../src)
extern "C" void SpecOrthSchur(double *MatR, double *VecA,
                              double *MatQ, int n);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // -------------------------------------------------------------
    // Check inputs / outputs
    // -------------------------------------------------------------
    if (nrhs != 1)
    {
        mexErrMsgIdAndTxt("skewblas:so2schur:nrhs",
                          "One input required: Q (special-orthogonal n x n).");
    }
    if (nlhs != 2)
    {
        mexErrMsgIdAndTxt("skewblas:so2schur:nlhs",
                          "Two outputs required: [R, A].");
    }

    const mxArray *Q_in = prhs[0];

    if (!mxIsDouble(Q_in) || mxIsComplex(Q_in))
    {
        mexErrMsgIdAndTxt("skewblas:so2schur:type",
                          "Input Q must be a real double matrix.");
    }

    mwSize n = mxGetM(Q_in);
    if (n != mxGetN(Q_in))
    {
        mexErrMsgIdAndTxt("skewblas:so2schur:square",
                          "Input Q must be square.");
    }

    // -------------------------------------------------------------
    // Allocate outputs
    // R : n x n real
    // A : (n/2) x 1 real
    // -------------------------------------------------------------
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n / 2, 1, mxREAL);

    double *R = mxGetPr(plhs[0]);
    double *A = mxGetPr(plhs[1]);
    double *Q = mxGetPr(Q_in);

    // -------------------------------------------------------------
    // Call backend C++ routine
    // -------------------------------------------------------------
    SpecOrthSchur(R, A, Q, static_cast<int>(n));
}
