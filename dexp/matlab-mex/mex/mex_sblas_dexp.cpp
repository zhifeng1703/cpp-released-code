#include "mex.h"
#include "matrix.h"
#include <cstring>

extern "C"
{

    void SkewSymmSchur(double *MatR, double *VecA, double *MatS, int n);

    // Compute both forward & inverse DK parameters
    void DexpSkewSymmPara(double *ParaForward, double *ParaInverse,
                          double *VecA, int n);

    void DexpSkewSymmForward(double *MatY, double *MatX,
                             double *MatR, double *ParaForward, int n);

    void DexpSkewSymmInverse(double *MatY, double *MatX,
                             double *MatR, double *ParaInverse, int n);
}

static void fail(const char *msg)
{
    mexErrMsgIdAndTxt("skewblas:mex_sblas_dexp", msg);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs < 3)
        fail("Usage: Y = mex_sblas_dexp(S, X, dir, [, para_struct])");

    /* ---------------------------------------------------------
       Read S
       --------------------------------------------------------- */
    const mxArray *S_mx = prhs[0];
    if (!mxIsDouble(S_mx) || mxIsComplex(S_mx))
        fail("S must be a real matrix.");

    int n = mxGetM(S_mx);
    if (mxGetN(S_mx) != n)
        fail("S must be n×n.");

    double *S = mxGetPr(S_mx);

    /* ---------------------------------------------------------
       Read X
       --------------------------------------------------------- */
    const mxArray *X_mx = prhs[1];
    if (!mxIsDouble(X_mx) || mxIsComplex(X_mx))
        fail("X must be a real matrix.");

    if (mxGetM(X_mx) != n || mxGetN(X_mx) != n)
        fail("X must be n×n.");

    double *X = mxGetPr(X_mx);

    /* ---------------------------------------------------------
       Read direction: 'fwd' or 'inv'
       --------------------------------------------------------- */
    char dirbuf[16];
    mxGetString(prhs[2], dirbuf, sizeof(dirbuf));
    bool forward = !(dirbuf[0] == 'i' || dirbuf[0] == 'I');

    /* ---------------------------------------------------------
       Check if user supplied para struct: para = struct(R,A,ParaForward,ParaInverse)
       --------------------------------------------------------- */
    bool have_para = (nrhs >= 4) && mxIsStruct(prhs[3]);

    mxArray *R_mx, *A_mx, *PF_mx, *PI_mx;
    double *R, *A, *ParaF, *ParaI;

    int a = n / 2;
    int par_size = a * (a - 1) * 8 + a * 4;

    if (!have_para)
    {
        /* ==============================================
           Compute Schur + DK parameters
           ============================================== */

        R_mx = mxCreateDoubleMatrix(n, n, mxREAL);
        A_mx = mxCreateDoubleMatrix(a, 1, mxREAL);

        R = mxGetPr(R_mx);
        A = mxGetPr(A_mx);

        SkewSymmSchur(R, A, S, n);

        PF_mx = mxCreateDoubleMatrix(par_size, 1, mxREAL);
        PI_mx = mxCreateDoubleMatrix(par_size, 1, mxREAL);

        ParaF = mxGetPr(PF_mx);
        ParaI = mxGetPr(PI_mx);

        DexpSkewSymmPara(ParaF, ParaI, A, n);
    }
    else
    {
        /* ==============================================
           Load provided para struct
           ============================================== */

        const mxArray *para_struct = prhs[3];

        R_mx = mxGetField(para_struct, 0, "R");
        A_mx = mxGetField(para_struct, 0, "A");
        PF_mx = mxGetField(para_struct, 0, "Fwd");
        PI_mx = mxGetField(para_struct, 0, "Inv");

        if (!R_mx || !A_mx || !PF_mx || !PI_mx)
            fail("para struct must contain fields R, A, Fwd, Inv");

        R = mxGetPr(R_mx);
        A = mxGetPr(A_mx);
        ParaF = mxGetPr(PF_mx);
        ParaI = mxGetPr(PI_mx);
    }

    /* ---------------------------------------------------------
       Output Y
       --------------------------------------------------------- */
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    double *Y = mxGetPr(plhs[0]);

    /* ---------------------------------------------------------
       Compute forward or inverse differential
       --------------------------------------------------------- */
    if (forward)
    {
        DexpSkewSymmForward(Y, X, R, ParaF, n);
    }
    else
    {
        DexpSkewSymmInverse(Y, X, R, ParaI, n);
    }
}
