// ===============================================================
//  skewblas_mex.cpp
//  MEX-only interface layer for MATLAB
// ===============================================================

#ifdef MATLAB_MEX_BUILD

#include "matOp.hpp"
#include "dexpSkew.hpp"
#include "dexpPade.hpp"
#include "dexpSkewEig.hpp"
#include "dlogSkewPade.hpp"
#include "skewSchFac.hpp"
#include "colMat.hpp"
#include <cstring>

// ===============================================================
//  ExpmPade
// ===============================================================
extern "C" void ExpmPade(double *MatQ, double *MatS, int n, int p)
{
    ColMat<double> TmpS(n, n);
    std::memcpy(TmpS.v, MatS, sizeof(double) * n * n);

    expmPadeApprox pade(n);
    pade.Expm(MatQ, n, TmpS.v, n, p);
}

// ===============================================================
//  ExpmSkewSymm
// ===============================================================
extern "C" void ExpmSkewSymm(double *MatQ, double *MatS, int n)
{
    SkewSchurFactor ssf(n);

    std::memcpy(ssf.H.MatH.v, MatS, sizeof(double) * n * n);
    ssf.SchurAngular_SkewSymm();
    ssf.Explict_Vector();
    ssf.Exponential(MatQ, n);
}

// ===============================================================
//  SkewSymmSchur
// ===============================================================
extern "C" void SkewSymmSchur(double *MatR, double *VecA, double *MatS, int n)
{
    ColMat<double> TmpS(n, n);
    SkewSchurFactor ssf(n);

    TmpS.Assign(MatS, n);

    std::memset(VecA, 0, sizeof(double) * (n / 2));

    auto saved_R = ssf.R.v;
    auto saved_A = ssf.A.v;

    ssf.R.v = MatR;
    ssf.A.v = VecA;
    ssf.a = n / 2;

    ssf.Factor_SkewSymm(TmpS.v, n);
    ssf.Explict_Vector();

    ssf.R.v = saved_R;
    ssf.A.v = saved_A;
}

// ===============================================================
//  SpecOrthSchur
// ===============================================================
extern "C" void SpecOrthSchur(double *MatR, double *VecA, double *MatQ, int n)
{
    ColMat<double> TmpQ(n, n);
    SkewSchurFactor ssf(n);

    TmpQ.Assign(MatQ, n);

    std::memset(VecA, 0, sizeof(double) * (n / 2));

    auto saved_R = ssf.R.v;
    auto saved_A = ssf.A.v;

    ssf.R.v = MatR;
    ssf.A.v = VecA;
    ssf.a = n / 2;

    ssf.Factor_SpecOrth(TmpQ.v, n);
    ssf.Explict_Vector();

    ssf.R.v = saved_R;
    ssf.A.v = saved_A;
}

// ===============================================================
//  BuildSpecOrthSchur
// ===============================================================
extern "C" void BuildSpecOrthSchur(double *MatQ, double *MatR, double *VecA, int n)
{
    SkewSchurFactor ssf(n);

    auto saved_R = ssf.R.v;
    auto saved_A = ssf.A.v;

    ssf.R.v = MatR;
    ssf.A.v = VecA;
    ssf.a = n / 2;

    ssf.Exponential(MatQ, n);

    ssf.R.v = saved_R;
    ssf.A.v = saved_A;
}

// ===============================================================
//  BuildSkewSymmSchur
// ===============================================================
extern "C" void BuildSkewSymmSchur(double *MatS, double *MatR, double *VecA, int n)
{
    SkewSchurFactor ssf(n);

    auto saved_R = ssf.R.v;
    auto saved_A = ssf.A.v;

    ssf.R.v = MatR;
    ssf.A.v = VecA;
    ssf.a = n / 2;

    ssf.GetSkewSymm(MatS, n);
    skewl2m(MatS, n, n);
    for (int i = 0; i < n; i += n + 1)
        MatS[i] = 0;

    ssf.R.v = saved_R;
    ssf.A.v = saved_A;
}

// ===============================================================
//  DexpSkewSymmPara
// ===============================================================
extern "C" void DexpSkewSymmPara(double *ParaForward, double *ParaInverse, double *VecA, int n)
{
    dexpSkewSymmPara dexp(n);

    auto saved_f = dexp._forward_para.v;
    auto saved_i = dexp._inverse_para.v;

    dexp._forward_para.v = ParaForward;
    dexp._inverse_para.v = ParaInverse;
    dexp.a = n / 2;
    dexp._set_A(VecA);

    dexp._forward_para.v = saved_f;
    dexp._inverse_para.v = saved_i;
}

// ===============================================================
//  DexpSkewSymmForward
// ===============================================================
extern "C" void DexpSkewSymmForward(double *MatY, double *MatX,
                                    double *MatR, double *ParaForward,
                                    int n)
{
    dexpSkewSymmPara dexp(n);
    SkewSymmMat SSX(n), SSY(n);

    auto saved_f = dexp._forward_para.v;
    auto saved_R = dexp.R.v;
    auto saved_vx = SSX.v;
    auto saved_vy = SSY.v;

    dexp._forward_para.v = ParaForward;
    dexp.R.v = MatR;
    SSX.v = MatX;
    SSY.v = MatY;
    SSX.init_low_vec();
    SSY.init_low_vec();
    SSX.mat2vec();

    dexp.Forward(SSY, SSX);

    dexp._forward_para.v = saved_f;
    dexp.R.v = saved_R;
    SSX.v = saved_vx;
    SSY.v = saved_vy;
}

// ===============================================================
//  DexpSkewSymmInverse
// ===============================================================
extern "C" void DexpSkewSymmInverse(double *MatY, double *MatX,
                                    double *MatR, double *ParaInverse,
                                    int n)
{
    dexpSkewSymmPara dexp(n);
    SkewSymmMat SSX(n), SSY(n);

    auto saved_i = dexp._inverse_para.v;
    auto saved_R = dexp.R.v;
    auto saved_vx = SSX.v;
    auto saved_vy = SSY.v;

    dexp._inverse_para.v = ParaInverse;
    dexp.R.v = MatR;
    SSX.v = MatX;
    SSY.v = MatY;
    SSX.init_low_vec();
    SSY.init_low_vec();
    SSX.mat2vec();

    dexp.Inverse(SSY, SSX);

    dexp._inverse_para.v = saved_i;
    dexp.R.v = saved_R;
    SSX.v = saved_vx;
    SSY.v = saved_vy;
}

// ===============================================================
//  Daletskii–Krein (Eigen-basis)
// ===============================================================
extern "C" void DexpDalKreinParaSkewSymm(double *CmpxParaForward,
                                         double *CmpxParaInverse,
                                         double *CmpxEigVal,
                                         int n)
{
    dexpSkewEigenPara dexp(n);

    auto saved_f = dexp._forward_para.v;
    auto saved_i = dexp._inverse_para.v;
    auto saved_e = dexp.EigVal.v;

    dexp._forward_para.v = reinterpret_cast<CMPX_TYPE *>(CmpxParaForward);
    dexp._inverse_para.v = reinterpret_cast<CMPX_TYPE *>(CmpxParaInverse);
    dexp.EigVal.v = reinterpret_cast<CMPX_TYPE *>(CmpxEigVal);

    dexp.Parameter_Skew();

    dexp._forward_para.v = saved_f;
    dexp._inverse_para.v = saved_i;
    dexp.EigVal.v = saved_e;
}

// ===============================================================
//  Daletskii–Krein full tangent
// ===============================================================
extern "C" void DexpDalKreinParaFulTange(double *CmpxParaForward,
                                         double *CmpxParaInverse,
                                         double *CmpxEigVal,
                                         int n)
{
    dexpSkewEigenPara dexp(n);

    auto saved_f = dexp._forward_para.v;
    auto saved_i = dexp._inverse_para.v;
    auto saved_e = dexp.EigVal.v;

    dexp._forward_para.v = reinterpret_cast<CMPX_TYPE *>(CmpxParaForward);
    dexp._inverse_para.v = reinterpret_cast<CMPX_TYPE *>(CmpxParaInverse);
    dexp.EigVal.v = reinterpret_cast<CMPX_TYPE *>(CmpxEigVal);

    dexp.Parameter_Full();

    dexp._forward_para.v = saved_f;
    dexp._inverse_para.v = saved_i;
    dexp.EigVal.v = saved_e;
}

// ===============================================================
//  DK forward differential
// ===============================================================
extern "C" void DexpDalKreinForward(double *RealMatY, double *RealMatX,
                                    double *CmpxMatEigVec,
                                    double *CmpxParaForward,
                                    int n)
{
    dexpSkewEigenPara dexp(n);

    auto saved_f = dexp._forward_para.v;
    auto saved_v = dexp.EigVec.v;

    dexp._forward_para.v =
        reinterpret_cast<CMPX_TYPE *>(CmpxParaForward);

    dexp.EigVec.v =
        reinterpret_cast<CMPX_TYPE *>(CmpxMatEigVec);

    dexp.Dexp(RealMatY, n, RealMatX, n);

    dexp._forward_para.v = saved_f;
    dexp.EigVec.v = saved_v;
}

// ===============================================================
//  DK inverse differential
// ===============================================================
extern "C" void DexpDalKreinInverse(double *RealMatY, double *RealMatX,
                                    double *CmpxMatEigVec,
                                    double *CmpxParaInverse,
                                    int n)
{
    dexpSkewEigenPara dexp(n);

    auto saved_i = dexp._inverse_para.v;
    auto saved_v = dexp.EigVec.v;

    dexp._inverse_para.v =
        reinterpret_cast<CMPX_TYPE *>(CmpxParaInverse);

    dexp.EigVec.v =
        reinterpret_cast<CMPX_TYPE *>(CmpxMatEigVec);

    dexp.DexpInv(RealMatY, n, RealMatX, n);

    dexp._inverse_para.v = saved_i;
    dexp.EigVec.v = saved_v;
}

// ===============================================================
//  DexpPade forward / inverse
// ===============================================================
extern "C" void DexpPadeForward(double *MatY, double *MatX,
                                double *MatA, int n, int p)
{
    dexpPadeApprox dexp(n);
    dexp.Expm(MatA, n, p);
    dexp.Dexp(MatY, n, MatX, n);
}

extern "C" void DexpPadeInverse(double *MatY, double *MatX,
                                double *MatR, double *VecA,
                                int dim, int pade_order, int scale_order)
{
    dlogSkewPadeApprox dlog(dim);

    dlog._set_R(MatR, dim);
    dlog._set_A(VecA);

    // int scale = std::ceil(std::log2(normS / _EXPM_PADE_APPROX_BOUNDS[12]));
    if (scale_order < 0)
        scale_order = 0;

    dlog._set_orders(pade_order, scale_order);

    dlog.Dlog(MatY, dim, MatX, dim);
}

#endif // MATLAB_MEX_BUILD
