#include "mex.h"
#include "matrix.h"
#include <string>
#include <cstring>

extern "C"
{
    void DexpDalKreinParaSkewSymm(double *CmpxParaForward, double *CmpxParaInverse,
                                  double *CmpxEigVal, int n);

    void DexpDalKreinParaFulTange(double *CmpxParaForward, double *CmpxParaInverse,
                                  double *CmpxEigVal, int n);
}

static void mexErr(const char *msg)
{
    mexErrMsgIdAndTxt("skewblas:dkpara", msg);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /*
        Usage:
            [PF, PI] = mex_sblas_dkpara(EigVal, mode)

        where:
            - EigVal : 2*n real vector (interleaved complex eigenvalues)
            - mode   : 'skew' or 'full'
    */

    if (nrhs != 2)
        mexErr("Usage: [PF, PI] = mex_sblas_dkpara(EigVal, mode)");

    /* ============================================================
       1) Read EigVal (2*n real interleaved vector)
       ============================================================ */
    const mxArray *EigVal_mx = prhs[0];
    if (!mxIsDouble(EigVal_mx) || mxIsComplex(EigVal_mx))
        mexErr("EigVal must be a real 2*n vector (interleaved).");

    mwSize lenEV = mxGetNumberOfElements(EigVal_mx);
    if (lenEV % 2 != 0)
        mexErr("EigVal length must be even (2*n).");

    int n = static_cast<int>(lenEV / 2);
    double *EV = mxGetPr(EigVal_mx);

    /* ============================================================
       2) Read mode keyword
       ============================================================ */
    char mode_str[32];
    mxGetString(prhs[1], mode_str, sizeof(mode_str));
    std::string mode(mode_str);

    bool use_full = false;
    if (mode == "skew")
        use_full = false;
    else if (mode == "full")
        use_full = true;
    else
        mexErr("mode must be 'skew' or 'full'.");

    /* ============================================================
       3) Allocate outputs PF, PI (2*n*n real)
       ============================================================ */
    plhs[0] = mxCreateDoubleMatrix(2 * n * n, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(2 * n * n, 1, mxREAL);

    double *PF = mxGetPr(plhs[0]);
    double *PI = mxGetPr(plhs[1]);

    /* ============================================================
       4) Call backend generator
       ============================================================ */
    if (use_full)
        DexpDalKreinParaFulTange(PF, PI, EV, n);
    else
        DexpDalKreinParaSkewSymm(PF, PI, EV, n);
}
