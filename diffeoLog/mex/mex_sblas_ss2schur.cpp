#include "mex.h"
#include "matrix.h"

// C interface provided by skewblas_mex.cpp (compiled from ../src)
extern "C" void SkewSymmSchur(double *MatR, double *VecA,
                              double *MatS, int n);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // -------------------------------------------------------------
    // Check inputs / outputs
    // -------------------------------------------------------------
    if (nrhs != 1)
    {
        mexErrMsgIdAndTxt("skewblas:ss2schur:nrhs",
                          "One input required: S (skew-symmetric n x n).");
    }
    if (nlhs != 2)
    {
        mexErrMsgIdAndTxt("skewblas:ss2schur:nlhs",
                          "Two outputs required: [R, A].");
    }

    const mxArray *S_in = prhs[0];

    if (!mxIsDouble(S_in) || mxIsComplex(S_in))
    {
        mexErrMsgIdAndTxt("skewblas:ss2schur:type",
                          "Input S must be a real double matrix.");
    }

    mwSize n = mxGetM(S_in);
    if (n != mxGetN(S_in))
    {
        mexErrMsgIdAndTxt("skewblas:ss2schur:square",
                          "Input S must be square.");
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
    double *S = mxGetPr(S_in); // read-only in our implementation

    // -------------------------------------------------------------
    // Call backend C++ routine
    // -------------------------------------------------------------
    SkewSymmSchur(R, A, S, static_cast<int>(n));
}
