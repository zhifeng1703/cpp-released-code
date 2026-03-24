#include "skewblas_api.hpp"

#include "dexpSkew.hpp"
#include "dexpPade.hpp"
#include "dexpSkewEig.hpp"
#include "dlogSkewPade.hpp"

extern "C"
{
    void SkewSymmSchur(double *MatR, double *VecA, double *MatS, int n)
    {
        auto ssf = SkewSchurFactor(n);

        ssf.Factor(MatS, n);
        ssf.Explict_Vector(MatR, n);
        memset(VecA, 0, sizeof(double) * (n / 2));
        memcpy(VecA, ssf.A, sizeof(double) * ssf.a);
    }

    void SpecOrthSchur(double *MatR, double *VecA, double *MatQ, int n)
    {
        auto ssf = SkewSchurFactor(n);

        // Evil implementation:

        memset(VecA, 0, sizeof(double) * (n / 2));
        auto temp_r = ssf.R.v, temp_a = ssf.A.v;

        ssf.R.v = MatR;
        ssf.A.v = VecA;

        ssf.Factor_SpecOrth(MatQ, n);

        ssf.R.v = temp_r;
        ssf.A.v = temp_a;

        // ssf.H.Skew_Assign(MatQ, n);
        // ssf.SchurAngular_SkewSymm();
        // ssf.Explict_Vector(MatR, n);
        // ssf._principal_angles(MatQ, n, MatR, n);
        // memset(VecA, 0, sizeof(double) * (n / 2));
        // memcpy(VecA, ssf.A, sizeof(double) * ssf.a);
    }

    // // ---------------------------------------------------------------------------
    // //  Forward differential exponential
    // // ---------------------------------------------------------------------------
    // void dexp_forward(double *A, double *M, double *N, int n)
    // {
    //     SkewSymmMat MatA(n), MatM(n), MatN(n), MatQ(n, n), MatY(n, n);
    //     MatA.v = A;
    //     MatM.v = M;
    //     MatN.v = N;

    //     auto ssf = SkewSchurFactor(n);
    //     auto dexp_skew = dexpSkewSymmPara(n);
    //     auto dexp_pade = dexpPadeApprox(n);
    //     auto dexp_cmpx = dexpSkewEigenPara(n);

    //     // 1. Schur factorization
    //     ssf.Factor(MatA.v, n);
    //     ssf.Explict_Vector();

    //     // 2. Skew-based exponential and directional exp
    //     dexp_skew.Parameter(ssf);
    //     dexp_skew.Forward(MatN, MatM);

    //     // 3. (Optionally) direct matrix exponential Q = exp(A)
    //     ssf.Exponential(MatQ.v, n);

    //     // 4. Padé approximation (order 13)
    //     dexp_pade.Expm(MatQ.v, n, MatA.v, n, 13);
    //     dexp_pade.Dexp(MatY.v, n, MatM.v, n);

    //     // 5. Semi-simple (eigen-based)
    //     dexp_cmpx.Assign(ssf);
    //     dexp_cmpx.Parameter();
    //     dexp_cmpx.Dexp(MatN.v, n, MatM.v, n);

    //     // Results written into N (and possibly A overwritten as exp(A))
    // }

    // // ---------------------------------------------------------------------------
    // //  Inverse differential logarithm
    // // ---------------------------------------------------------------------------
    // void dlog_inverse(double *A, double *M, double *N, int n)
    // {
    //     SkewSymmMat MatA(n), MatM(n), MatN(n), MatQ(n, n), MatY(n, n);
    //     MatA.v = A;
    //     MatM.v = M;
    //     MatN.v = N;

    //     auto ssf = SkewSchurFactor(n);
    //     auto dlog_pade = dlogSkewPadeApprox(n);
    //     auto dexp_cmpx = dexpSkewEigenPara(n);

    //     // Factorize and prepare
    //     ssf.Factor(MatA.v, n);
    //     ssf.Explict_Vector();
    //     dlog_pade.Parameter(ssf);

    //     // Determine scaling for stability
    //     int scale = ceil(log2(MatA.Norm1() / _EXPM_PADE_APPROX_BOUNDS[12]));
    //     scale = std::max(scale, 0);
    //     dlog_pade._set_orders(7, scale);

    //     // Compute exponential once
    //     ssf.Exponential(MatQ.v, n);

    //     // Combine
    //     cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    //                 n, n, n, 1.0, MatQ.v, n, MatN.v, n, 0.0, MatY.v, n);

    //     // Apply differential log
    //     dlog_pade.Dlog(MatN.v, n, MatY.v, n);

    //     // Semisimple inverse
    //     dexp_cmpx.Assign(ssf);
    //     dexp_cmpx.Parameter();
    //     dexp_cmpx.DexpInv(MatN.v, n, MatM.v, n);
    // }

} // extern "C"
