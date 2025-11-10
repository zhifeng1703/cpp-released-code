#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

    void SkewSymmSchur(double *MatR, double *VecA, double *MatS, int n);

    void SpecOrthSchur(double *MatR, double *VecA, double *MatQ, int n);

    // void SkewSymmDexpPara(double *)

    // /**
    //  * Forward directional exponential on a skew-symmetric matrix.
    //  *   Computes Q = exp(A), N = d exp(A)[M].
    //  * Arguments:
    //  *   A  – input skew-symmetric matrix (n×n, column-major)
    //  *   M  – perturbation matrix (n×n)
    //  *   N  – output matrix (n×n)
    //  *   n  – dimension
    //  */
    // void dexp_forward(double *A, double *M, double *N, int n);

    // /**
    //  * Inverse directional logarithm on a skew-symmetric matrix.
    //  *   Computes N = d log(A)[M].
    //  */
    // void dlog_inverse(double *A, double *M, double *N, int n);

#ifdef __cplusplus
}
#endif
