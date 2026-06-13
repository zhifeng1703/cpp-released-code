#include "mex.h"
#include "matrix.h"
#include <cstring>

extern "C"
{
    void DexpSkewSymmPara(double *ParaForward, double *ParaInverse,
                          double *A, int m);
}

static void fail(const char *msg)
{
    mexErrMsgIdAndTxt("skewblas:mex_sblas_dexp_para", msg);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* ------------------------------------------------------
       Check input
       ------------------------------------------------------ */
    if (nrhs != 1)
        fail("Usage: [ParaF, ParaI] = mex_sblas_dexp_para(A)");

    const mxArray *A_mx = prhs[0];

    if (!mxIsDouble(A_mx) || mxIsComplex(A_mx))
        fail("Input A must be a real double vector.");

    int m = (int)mxGetNumberOfElements(A_mx);
    if (m <= 0)
        fail("Input A must be nonempty.");

    double *A = mxGetPr(A_mx);

    /* ------------------------------------------------------
       Prepare output parameter vectors
       ------------------------------------------------------ */

    int par_size = m * (m - 1) * 8 + m * 4;
    plhs[0] = mxCreateDoubleMatrix(par_size, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(par_size, 1, mxREAL);

    double *ParaF = mxGetPr(plhs[0]);
    double *ParaI = mxGetPr(plhs[1]);

    /* ------------------------------------------------------
       Compute the forward & inverse DK parameters
       ------------------------------------------------------ */
    DexpSkewSymmPara(ParaF, ParaI, A, 2 * m);
}
