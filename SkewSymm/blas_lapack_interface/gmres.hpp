//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the
// SIAM Templates book.
//
// The return value indicates convergence within max_iterTotal(input)
// iterations (0), or no convergence within max_iterTotaliterations (1).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iterTotal --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//*****************************************************************

// Zhifeng Deng: This implementation is modified from the template provided by netlib
// https://www.netlib.org/templates/cpp//gmres.h

#pragma once

#include <cmath>

#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"
#include "matOp.hpp"
#include "vecOp.hpp"

#include "ddot.hpp"
#include "daxpy.hpp"

// Preconditioner version not implemented.

REAL_TYPE *_GMRES_Work(INTE_TYPE n, INTE_TYPE m);

REAL_TYPE _GMRES_NormEst(REAL_TYPE alpha, INTE_TYPE k, ColMat<REAL_TYPE> &H, REAL_TYPE *s);

void _GMRES_Update(REAL_TYPE *x, INTE_TYPE n, INTE_TYPE k, ColMat<REAL_TYPE> &H, REAL_TYPE *s, ColMat<REAL_TYPE> &V, REAL_TYPE *w);

void _GMRES_Generate_Rotation(REAL_TYPE &dx, REAL_TYPE &dy, REAL_TYPE &cs, REAL_TYPE &sn);

void _GMRES_Rotate(REAL_TYPE &dx, REAL_TYPE &dy, REAL_TYPE &cs, REAL_TYPE &sn);

template <class Operator>
INTE_TYPE GMRES(const Operator &A, REAL_TYPE *x, const REAL_TYPE *b, ColMat<REAL_TYPE> &H,
                INTE_TYPE n, INTE_TYPE &MaxIter, INTE_TYPE &m, REAL_TYPE &RelTol, REAL_TYPE &MaxNorm,
                REAL_TYPE *workspace, ColMat<REAL_TYPE> &V)
{
    // This implementation solves the matrix-free linear system Ax = b, where A is n x n.
    // At the k-th iteration, a solution x_k in the Krylov space
    // Krylov(k) := col(r, A * r, A^2 * r, ..., A^k * r) where r = Ax_0 - b
    // is found such that it minimizes x_k = \argmin_{x\in Krylov(k)} |Ax - b|.
    // At the n-th iteration, the solution to Ax = b is guaranteed.

    // At the k-th iteration, a call to the action Ax and an extra O(kn) FLOPs are needed.
    // The parameter m control to growth of k to avoid growing complexity.
    // At the m-th iteration, if the current attempt x_m does not meet the stop criteria
    // GMRES restarts by settting x_0 <- x_m, r <- b - A * x_0, k <- 0.
    // The parameter MaxIter controls the overall number of iterations, which is also
    // the total number of all to the action A * x.

    // GMRES requires workspace of H[m + 1, m] for the Upper Hessenberg matrix, V[n, m + 1] for the Krylov basis,
    // and w[2n + 3m +3] = (r[n], y[n], s[m+1], cs[m+1], sn[m+1]) for temporary vectors.
    REAL_TYPE resid;
    INTE_TYPE iter, iterTotal = 1, krylov;
    // workspace = (r[n], y[n], s[m+1], cs[m+1], sn[m+1])
    REAL_TYPE *r = workspace;
    REAL_TYPE *y = r + n;
    REAL_TYPE *s = y + n;
    REAL_TYPE *cs = s + m + 1;
    REAL_TYPE *sn = cs + m + 1;

    REAL_TYPE normb = norml2(n, b);
    memcpy(r, b, sizeof(REAL_TYPE) * n);
    A.Remainder(r, x);
    REAL_TYPE alpha = norml2(n, x);
    REAL_TYPE beta = norml2(n, r);

    if (normb == 0.0)
        normb = 1;

    if ((resid = norml2(n, r) / normb) <= RelTol)
    {
        RelTol = resid;
        MaxIter = 0;
        return 0;
    }

    // Vector *v = new Vector[m+1];

    // ColMat<REAL_TYPE> V = ColMat<REAL_TYPE>(n, m+1);
    // V.fast_col_access();

    while (iterTotal <= MaxIter)
    {
        // v[0] = r * (1.0 / beta);                                 // V[:, 0] = r / beta
        memcpy(V.fcol(0), r, sizeof(REAL_TYPE) * n);
        scal(n, 1.0 / beta, V.fcol(0));

        memset(s, 0, sizeof(REAL_TYPE) * (m + 1));
        s[0] = beta;

        // printf("\tGMRES Outer Loop Iteration\t Vector norm of current x:\t %1.8f, \t Vector norm of current residual:\t %1.8f\n", alpha, beta);

        for (iter = 0; iter < m && iterTotal <= MaxIter; iter++, iterTotal++)
        {
            // printf("\t\tGMRES Inner Loop Iteration:\t %lld,\t Vector norm of current residual:\t %1.8f\n", iterTotal, abs(s[iter]));

            A.Action(y, V.fcol(iter)); // y = A * V[:, i]
            for (krylov = 0; krylov <= iter; krylov++)
            {
                H.fvar(krylov, iter) = my_ddot(n, y, V.fcol(krylov));
                my_daxpy(n, -H.fvar(krylov, iter), V.fcol(krylov), y); // y = y - H[k, i] * V[:, k], k = 0, ..., i
            }
            H.fvar(iter + 1, iter) = norml2(n, y);                       // H[i+1, i] = norm(y)
            memcpy(V.fcol(iter + 1), y, sizeof(REAL_TYPE) * n);          // V[:, iter+ 1] = y / H[i+1, i]
            my_dscal(n, 1.0 / H.fvar(iter + 1, iter), V.fcol(iter + 1)); // V[:, iter+ 1] = y / H[i+1, i]

            for (krylov = 0; krylov < iter; krylov++)
                _GMRES_Rotate(H.fvar(krylov, iter), H.fvar(krylov + 1, iter), cs[krylov], sn[krylov]);

            _GMRES_Generate_Rotation(H.fvar(iter, iter), H.fvar(iter + 1, iter), cs[iter], sn[iter]);
            _GMRES_Rotate(H.fvar(iter, iter), H.fvar(iter + 1, iter), cs[iter], sn[iter]);
            _GMRES_Rotate(s[iter], s[iter + 1], cs[iter], sn[iter]);

            if ((resid = abs(s[iter + 1]) / normb) < RelTol)
            {
                _GMRES_Update(x, n, iter, H, s, V, y);
                RelTol = resid;
                MaxIter = iterTotal;
                // printf("\tGMRES exit at iteration:\t %lld,\t with relative error:\t %1.8f\n", iterTotal, resid);

                return 0;
            }
        }
        _GMRES_Update(x, n, m - 1, H, s, V, y);

        memcpy(r, b, sizeof(REAL_TYPE) * n);
        A.Remainder(r, x); // r = b - Ax

        alpha = norml2(n, x);
        if (MaxNorm != 0.0 && alpha > MaxNorm)
        {
            RelTol = resid;
            MaxIter = iterTotal;
            // printf("\tGMRES terminated at iteration:\t %lld,\t with huge norm of vector x:\t %1.8f\n", iterTotal, alpha);
            return 1;
        }

        beta = norml2(n, r);
        if ((resid = beta / normb) < RelTol)
        {
            RelTol = resid;
            MaxIter = iterTotal;
            // printf("\tGMRES exit at iteration:\t %lld,\t with relative error:\t %1.8f\n", iterTotal, resid);
            return 0;
        }
    }

    RelTol = resid;
    // printf("\tGMRES terminated at the maximum iteration\t %lld,\t with relative error:\t %1.8f\n", iterTotal, resid);
    return 1;
};
