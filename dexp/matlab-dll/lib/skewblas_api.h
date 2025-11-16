#ifndef SKEWBLAS_API_H
#define SKEWBLAS_API_H

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef _WIN32
#ifdef BUILDING_SKEWBLAS_DLL
#define SKEWBLAS_API __declspec(dllexport)
#else
#define SKEWBLAS_API __declspec(dllimport)
#endif
#else
#define SKEWBLAS_API
#endif

   /* ============================================================
      API functions for MATLAB / C interfaces
      ============================================================ */

   /**
    * ExpmPade
    *  Compute matrix exponential using Pade approximation.
    *
    *  Parameters:
    *   - MatQ [out] : n×n output matrix (double[n*n])
    *   - VecA [in]  : input flattened matrix (double[n*n])
    *   - n          : matrix dimension
    *   - p          : Pade order
    */
   SKEWBLAS_API void ExpmPade(double *MatQ,
                              double *VecA,
                              int n,
                              int p);

   /**
    * ExpmSkewSymm
    *  Compute exponential of skew-symmetric matrix via Schur decomposition.
    *
    *  Parameters:
    *   - MatQ [out] : n×n orthogonal exponential matrix
    *   - MatR [out] : n×n orthogonal matrix from factorization
    *   - VecA [out] : n/2 angular vector
    *   - MatS [in]  : n×n skew-symmetric matrix
    *   - n          : matrix dimension
    */
   SKEWBLAS_API void ExpmSkewSymm(double *MatQ,
                                  double *MatS,
                                  int n);

   /**
    * SkewSymmSchur
    *  Factorize skew-symmetric matrix MatS into MatR and VecA.
    */
   SKEWBLAS_API void SkewSymmSchur(double *MatR,
                                   double *VecA,
                                   double *MatS,
                                   int n);

   /**
    * SpecOrthSchur
    *  Factorize a special orthogonal matrix MatQ into MatR and VecA.
    */
   SKEWBLAS_API void SpecOrthSchur(double *MatR,
                                   double *VecA,
                                   double *MatQ,
                                   int n);

   /**
    * BuildSpecOrthSchur
    *  Reconstruct special-orthogonal matrix MatQ from MatR and VecA.
    */
   SKEWBLAS_API void BuildSpecOrthSchur(double *MatQ,
                                        double *MatR,
                                        double *VecA,
                                        int n);

   /**
    * BuildSkewSymmSchur
    *  Reconstruct skew-symmetric matrix MatS from MatR and VecA.
    */
   SKEWBLAS_API void BuildSkewSymmSchur(double *MatS,
                                        double *MatR,
                                        double *VecA,
                                        int n);

   /**
    * DexpSkewSymmPara
    *  Precompute forward/inverse parameters for skew-symmetric differential exponential.
    */
   SKEWBLAS_API void DexpSkewSymmPara(double *ParaForward,
                                      double *ParaInverse,
                                      double *VecA,
                                      int n);

   /**
    * DexpSkewSymmForward
    *  Apply differential exponential (forward) for skew-symmetric structure.
    */
   SKEWBLAS_API void DexpSkewSymmForward(double *MatY,
                                         double *MatX,
                                         double *MatR,
                                         double *ParaForward,
                                         int n);

   /**
    * DexpSkewSymmInverse
    *  Apply inverse differential exponential for skew-symmetric structure.
    */
   SKEWBLAS_API void DexpSkewSymmInverse(double *MatY,
                                         double *MatX,
                                         double *MatR,
                                         double *ParaInverse,
                                         int n);

   /**
    * DexpDalKreinParaSkewSymm
    *  Compute Daletskii–Kreĭn parameters for skew-symmetric case.
    */
   SKEWBLAS_API void DexpDalKreinParaSkewSymm(double *CmpxParaForward,
                                              double *CmpxParaInverse,
                                              double *CmpxEigVal,
                                              int n);

   /**
    * DexpDalKreinParaFulTange
    *  Compute Daletskii–Kreĭn parameters for full tangent case.
    */
   SKEWBLAS_API void DexpDalKreinParaFulTange(double *CmpxParaForward,
                                              double *CmpxParaInverse,
                                              double *CmpxEigVal,
                                              int n);

   /**
    * DexpDalKreinForward
    *  Compute forward Daletskii–Kreĭn differential.
    */
   SKEWBLAS_API void DexpDalKreinForward(double *RealMatY,
                                         double *RealMatX,
                                         double *CmpxMatEigVec,
                                         double *CmpxParaForward,
                                         int n);

   /**
    * DexpDalKreinInverse
    *  Compute inverse Daletskii–Kreĭn differential.
    */
   SKEWBLAS_API void DexpDalKreinInverse(double *RealMatY,
                                         double *RealMatX,
                                         double *CmpxMatEigVec,
                                         double *CmpxParaInverse,
                                         int n);

   /**
    * DexpPadeForward
    *  Compute forward differential of exponential via Pade approximation.
    */
   SKEWBLAS_API void DexpPadeForward(double *MatY,
                                     double *MatX,
                                     double *MatA,
                                     int n,
                                     int p);

   /**
    * DexpPadeInverse
    *  Compute inverse differential (logarithmic) via Pade approximation.
    */
   SKEWBLAS_API void DexpPadeInverse(double *MatY,
                                     double *MatX,
                                     double *MatR,
                                     double *VecA,
                                     double normS,
                                     int n,
                                     int p);

#ifdef __cplusplus
}
#endif

#endif /* SKEWBLAS_API_H */
