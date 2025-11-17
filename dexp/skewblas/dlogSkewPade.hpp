#pragma once

#include <fstream>

#include "colMat.hpp"
#include "skewSchFac.hpp"
#include "pivotedLU.hpp"

// This code implements the scaling and squaring method of the matrix exponential proposed in [Higham09].
// Not that it is the olde version that does not contained special treatment of the scaling order.
// For the improved version of the scaling and squaring mehtod, see expmss_improved.hpp.

// [Higham09]: Higham, Nicholas J. "The scaling and squaring method for the matrix exponential revisited."
// SIAM review 51, no. 4 (2009): 747-764.

// Note that this code serve as the preprocessing step (in the worst case scenario) of the directional
// derivative dexp implemented in dexpPadeSeries, which requires [13/13] Pade approximant.
// Therefore, this code always assumes the [13/13] Pade approximant is in use, without preprocessing on
// reducing the matrix 1-norm, i.e., it implements line 17 - 21 in [Higham09](Algorithm 2.3)

// #define _EXPM_PADE_APPROX_MAX_SCALING 10
// #define _EXPM_PADE_APPROX_MAX_PARANUM 13

// const INTE_TYPE _EXPM_PADE_APPROX_3[4] = {120, 60, 12, 1};
// const INTE_TYPE _EXPM_PADE_APPROX_5[6] = {30240, 15120, 3360, 420, 30, 1};
// const INTE_TYPE _EXPM_PADE_APPROX_7[8] = {17297280, 8652600, 1995840, 277200, 25200, 1512, 56, 1};
// const INTE_TYPE _EXPM_PADE_APPROX_9[10] = {17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1};
// const INTE_TYPE _EXPM_PADE_APPROX_13[14] = {64764752532480000, 32382376266240000, 7771770303897600,
//                                             1187353796428800, 129060195264000, 10559470521600,
//                                             670442572800, 33522128640, 1323241920,
//                                             40840800, 960960, 16380, 182, 1};

// const REAL_TYPE _EXPM_PADE_APPROX_BOUNDS[13] = {2.11e-8, 3.56e-4, 1.08e-2, 6.49e-2, 2.00e-1, 4.37e-1, 7.83e-1, 1.23e0, 1.78e0, 2.42e0, 3.13e0, 3.90e0, 4.74e0};

const REAL_TYPE _LOGM_PADE_ALPHA_7[] = {0.20897959, 0.19091503, 0.19091503, 0.13985270, 0.13985270, 0.06474248, 0.06474248};
const REAL_TYPE _LOGM_PADE_ALPHA_6[] = {0.23395697, 0.23395697, 0.18038079, 0.18038079, 0.08566225, 0.08566225};
const REAL_TYPE _LOGM_PADE_ALPHA_5[] = {0.28444444, 0.23931434, 0.23931434, 0.11846344, 0.11846344};
const REAL_TYPE _LOGM_PADE_ALPHA_4[] = {0.17392742, 0.17392742, 0.32607258, 0.32607258};
const REAL_TYPE _LOGM_PADE_ALPHA_3[] = {0.44444444, 0.27777778, 0.27777778};
const REAL_TYPE _LOGM_PADE_ALPHA_2[] = {0.50000000, 0.50000000};
const REAL_TYPE _LOGM_PADE_ALPHA_1[] = {1.00000000};
const REAL_TYPE *_LOGM_PADE_ALPHA[] = {_LOGM_PADE_ALPHA_1, _LOGM_PADE_ALPHA_2, _LOGM_PADE_ALPHA_3, _LOGM_PADE_ALPHA_4, _LOGM_PADE_ALPHA_5, _LOGM_PADE_ALPHA_6, _LOGM_PADE_ALPHA_7};

const REAL_TYPE _LOGM_PADE_BETA_7[] = {0.50000000, 0.70292258, 0.29707742, 0.87076559, 0.12923481, 0.97455396, 0.02544604};
const REAL_TYPE _LOGM_PADE_BETA_6[] = {0.61930959, 0.38069041, 0.83060469, 0.16939531, 0.96623476, 0.03376524};
const REAL_TYPE _LOGM_PADE_BETA_5[] = {0.50000000, 0.76923466, 0.23076534, 0.95308992, 0.04691008};
const REAL_TYPE _LOGM_PADE_BETA_4[] = {0.93056816, 0.06943184, 0.66999052, 0.33000948};
const REAL_TYPE _LOGM_PADE_BETA_3[] = {0.50000000, 0.88729833, 0.11270167};
const REAL_TYPE _LOGM_PADE_BETA_2[] = {0.78867513, 0.21132487};
const REAL_TYPE _LOGM_PADE_BETA_1[] = {0.50000000};
const REAL_TYPE *_LOGM_PADE_BETA[] = {_LOGM_PADE_BETA_1, _LOGM_PADE_BETA_2, _LOGM_PADE_BETA_3, _LOGM_PADE_BETA_4, _LOGM_PADE_BETA_5, _LOGM_PADE_BETA_6, _LOGM_PADE_BETA_7};

const REAL_TYPE _LOGM_PADE_APPROX_BOUNDS[] = {2.11e-8, 3.56e-4, 1.08e-2, 6.49e-2, 2.00e-1, 4.37e-1, 7.83e-1, 1.23e0, 1.78e0, 2.42e0, 3.13e0, 3.90e0, 4.74e0};

class dlogSkewPadeApprox
{
    typedef dlogSkewPadeApprox SELF_TYPE;
    typedef ArrVec<REAL_TYPE> VECT_TYPE;
    typedef ColMat<REAL_TYPE> MATX_TYPE;

public:
    INTE_TYPE d;
    INTE_TYPE s;
    INTE_TYPE m;

    const REAL_TYPE *alpha;
    const REAL_TYPE *beta;
    // Pre-computed constant type coefficients, not managed by this object.

    VECT_TYPE A;
    MATX_TYPE R;

    MATX_TYPE X;
    MATX_TYPE Y;

    MATX_TYPE Work;

    dlogSkewPadeApprox() : s(0), d(0), m(0), alpha(nullptr), beta(nullptr), A(VECT_TYPE()), R(MATX_TYPE()), X(MATX_TYPE()), Y(MATX_TYPE()), Work(MATX_TYPE()) {};
    dlogSkewPadeApprox(INTE_TYPE dim) : s(0), d(dim), m(0),
                                        A(VECT_TYPE(dim / 2)), R(MATX_TYPE(dim, dim)), X(MATX_TYPE(dim, dim)), Y(MATX_TYPE(dim, dim)), Work(MATX_TYPE(dim, dim + 1)) {};
    void copy(const SELF_TYPE &src)
    {
        this->alpha = src.alpha;
        this->beta = src.beta;
        this->A.copy(src.A);
        this->R.copy(src.R);
        this->X.copy(src.X);
        this->Y.copy(src.Y);
    }
    dlogSkewPadeApprox(const SELF_TYPE &src) : s(src.s), d(src.d), m(src.m),
                                               A(VECT_TYPE(src.d / 2)), R(MATX_TYPE(src.d, src.d)), X(MATX_TYPE(src.d, src.d)), Y(MATX_TYPE(src.d, src.d)),
                                               Work(MATX_TYPE(src.d, src.d + 1)) { this->copy(src); };
    void swap(SELF_TYPE &src)
    {
        this->A.swap(src.A);
        this->R.swap(src.R);
        this->X.swap(src.X);
        this->Y.swap(src.Y);
        this->Work.swap(src.Work);
        using std::swap;
        swap(this->d, src.d);
        swap(this->s, src.s);
        swap(this->m, src.m);
        swap(this->alpha, src.alpha);
        swap(this->beta, src.beta);
    }
    SELF_TYPE &operator=(const SELF_TYPE &rhs)
    {
        SELF_TYPE temp = SELF_TYPE(rhs);
        swap(temp);
        return (*this);
    }
    ~dlogSkewPadeApprox() {};

    void printMat(const char *name, const double *A, int row, int col)
    {
        FILE *f = fopen("debug.txt", "a"); // append mode
        if (!f)
            return;

        fprintf(f, "%s = [\n", name);

        for (int i = 0; i < row; ++i) // row index
        {
            for (int j = 0; j < col; ++j) // col index
            {
                fprintf(f, "% .15g ", A[i + j * row]); // column-major indexing
            }
            fprintf(f, ";\n");
        }

        fprintf(f, "]\n\n");
        fclose(f);
    };

    void printVec(const char *name, const double *A, int n)
    {
        FILE *f = fopen("debug.txt", "a"); // append mode
        if (!f)
            return;

        fprintf(f, "%s = [\n", name);

        for (int j = 0; j < n; ++j) // col index
        {
            fprintf(f, "% .15g ", A[j]); // column-major indexing
        }
        fprintf(f, ";\n");

        fprintf(f, "]\n\n");
        fclose(f);
    };

    void printMsg(const char *str)
    {
        FILE *f = fopen("debug.txt", "a"); // append mode
        if (!f)
            return;

        fprintf(f, "Dlog Setup: dim %d, pade order %d, scale order %d, %s.\n", d, m, s, str);
        fclose(f);
    };

    void Assign(REAL_TYPE *MatA, INTE_TYPE lda) { X.Assign(MatA, lda); };
    void Assign(const View_ColMat<REAL_TYPE> &ViewA) { Assign(ViewA.v, ViewA.ld); };
    void Assign(const ColMat<REAL_TYPE> &MatA) { Assign(MatA.v, MatA.r); };

    void _set_orders(INTE_TYPE order, INTE_TYPE scale)
    {
        s = scale;
        m = order;
        alpha = _LOGM_PADE_ALPHA[m - 1];
        beta = _LOGM_PADE_BETA[m - 1];
    };

    // Lyapunov Solver:

    // Rotation transpose R(-phi)
    static inline void rotT(REAL_TYPE Rm[4], REAL_TYPE phi)
    {
        REAL_TYPE c = std::cos(phi);
        REAL_TYPE s = std::sin(phi);
        Rm[0] = c;
        Rm[1] = s;
        Rm[2] = -s;
        Rm[3] = c;
    }

    // (R(phi)+I)^{-1} * b  (2×1)
    static inline void solve_rot_plus_I_left(REAL_TYPE y[2],
                                             const REAL_TYPE b[2],
                                             REAL_TYPE phi)
    {
        REAL_TYPE c = std::cos(phi);
        REAL_TYPE s = std::sin(phi);
        REAL_TYPE den = 2 * (1 + c);

        y[0] = ((1 + c) * b[0] + s * b[1]) / den;
        y[1] = (-(s)*b[0] + (1 + c) * b[1]) / den;
    }

    // b * (R(phi)+I)^{-1}  (1×2)
    static inline void solve_rot_plus_I_right(REAL_TYPE y[2],
                                              const REAL_TYPE b[2],
                                              REAL_TYPE phi)
    {
        REAL_TYPE c = std::cos(phi);
        REAL_TYPE s = std::sin(phi);
        REAL_TYPE den = 2 * (1 + c);

        y[0] = ((1 + c) * b[0] - s * b[1]) / den;
        y[1] = (s * b[0] + (1 + c) * b[1]) / den;
    }

    // 2×2 rotation–rotation symbolic formula:
    // Solve R(phi_p) X + X R(phi_q) = B
    static inline void solve_rot_rot_2x2(REAL_TYPE X[4],
                                         const REAL_TYPE B[4],
                                         REAL_TYPE phi_p,
                                         REAL_TYPE phi_q)
    {
        REAL_TYPE cp = std::cos(phi_p), sp = std::sin(phi_p);
        REAL_TYPE cq = std::cos(phi_q), sq = std::sin(phi_q);

        REAL_TYPE Delta = cp + cq;
        REAL_TYPE den = Delta * Delta + sp * sp + sq * sq; // = 2 + 2 cp cq

        REAL_TYPE RTp[4], RTq[4];
        rotT(RTp, phi_p);
        rotT(RTq, phi_q);

        REAL_TYPE BRtq[4], RtpB[4];

        // B * R(-phi_q)
        BRtq[0] = B[0] * RTq[0] + B[1] * RTq[2];
        BRtq[1] = B[0] * RTq[1] + B[1] * RTq[3];
        BRtq[2] = B[2] * RTq[0] + B[3] * RTq[2];
        BRtq[3] = B[2] * RTq[1] + B[3] * RTq[3];

        // R(-phi_p) * B
        RtpB[0] = RTp[0] * B[0] + RTp[1] * B[2];
        RtpB[1] = RTp[0] * B[1] + RTp[1] * B[3];
        RtpB[2] = RTp[2] * B[0] + RTp[3] * B[2];
        RtpB[3] = RTp[2] * B[1] + RTp[3] * B[3];

        // Combine symbolic terms
        X[0] = (Delta * B[0] + sp * BRtq[0] + sq * RtpB[0]) / den;
        X[1] = (Delta * B[1] + sp * BRtq[1] + sq * RtpB[1]) / den;
        X[2] = (Delta * B[2] + sp * BRtq[2] + sq * RtpB[2]) / den;
        X[3] = (Delta * B[3] + sp * BRtq[3] + sq * RtpB[3]) / den;
    }

    // MAIN BLOCK-LYAP FUNCTION
    // Solves A Y + Y A = X
    // where A = diag(R(theta[0]),...,R(theta[m-1]), 1?) with n=2m or 2m+1.
    void blk_lyap(REAL_TYPE *Y,
                  const REAL_TYPE *X,
                  const REAL_TYPE *theta)
    {
        INTE_TYPE n = d, nblocks = n / 2;
        bool has_scalar = (n == 2 * nblocks + 1);
        INTE_TYPE sidx = has_scalar ? (2 * nblocks) : -1;

        // Zero Y
        for (INTE_TYPE i = 0; i < n * n; i++)
            Y[i] = 0.0;

        // ================================
        //  1) Rotation–Rotation blocks (2×2)
        // ================================
        for (INTE_TYPE p = 0; p < nblocks; p++)
        {
            INTE_TYPE rp0 = 2 * p;
            INTE_TYPE rp1 = rp0 + 1;

            for (INTE_TYPE q = 0; q < nblocks; q++)
            {
                INTE_TYPE cq0 = 2 * q;
                INTE_TYPE cq1 = cq0 + 1;

                REAL_TYPE B[4];
                B[0] = X[rp0 + cq0 * n];
                B[1] = X[rp0 + cq1 * n];
                B[2] = X[rp1 + cq0 * n];
                B[3] = X[rp1 + cq1 * n];

                REAL_TYPE Z[4];
                solve_rot_rot_2x2(Z, B, theta[p], theta[q]);

                Y[rp0 + cq0 * n] = Z[0];
                Y[rp0 + cq1 * n] = Z[1];
                Y[rp1 + cq0 * n] = Z[2];
                Y[rp1 + cq1 * n] = Z[3];
            }
        }

        // ================================
        // Scalar block exists only if n=2m+1
        // ================================
        if (has_scalar)
        {
            // ------------------------------
            // 2) Rotation–Scalar (2×1)
            // ------------------------------
            for (INTE_TYPE p = 0; p < nblocks; p++)
            {
                INTE_TYPE rp0 = 2 * p;
                INTE_TYPE rp1 = rp0 + 1;

                REAL_TYPE b[2];
                b[0] = X[rp0 + sidx * n];
                b[1] = X[rp1 + sidx * n];

                REAL_TYPE yv[2];
                solve_rot_plus_I_left(yv, b, theta[p]);

                Y[rp0 + sidx * n] = yv[0];
                Y[rp1 + sidx * n] = yv[1];
            }

            // ------------------------------
            // 3) Scalar–Rotation (1×2)
            // ------------------------------
            for (INTE_TYPE q = 0; q < nblocks; q++)
            {
                INTE_TYPE cq0 = 2 * q;
                INTE_TYPE cq1 = cq0 + 1;

                REAL_TYPE b[2];
                b[0] = X[sidx + cq0 * n];
                b[1] = X[sidx + cq1 * n];

                REAL_TYPE yv[2];
                solve_rot_plus_I_right(yv, b, theta[q]);

                Y[sidx + cq0 * n] = yv[0];
                Y[sidx + cq1 * n] = yv[1];
            }

            // ------------------------------
            // 4) Scalar–Scalar (1×1)
            // ------------------------------
            REAL_TYPE bss = X[sidx + sidx * n];
            Y[sidx + sidx * n] = 0.5 * bss;
        }
    }

    // Inverse action

    // ------------------------------------------------------------
    // Compute 2x2 F = (I + beta*(R(phi)-I))^{-1} in closed form.
    // R(phi) = [c -s; s c],  K = R - I,  I + beta*K = [a -b; b a]
    // with a = 1 + beta(c-1), b = beta*s
    // => F = 1/(a^2 + b^2) [ a  b; -b  a ]
    // ------------------------------------------------------------
    static inline void compute_F_block(REAL_TYPE F[4],
                                       REAL_TYPE phi,
                                       REAL_TYPE beta)
    {
        REAL_TYPE c = std::cos(phi);
        REAL_TYPE s = std::sin(phi);

        REAL_TYPE a = 1.0 + beta * (c - 1.0);
        REAL_TYPE b = beta * s;
        REAL_TYPE den = a * a + b * b;

        REAL_TYPE inv = 1.0 / den;

        F[0] = a * inv;  // (1,1)
        F[1] = -b * inv; // (1,2)
        F[2] = b * inv;  // (2,1)
        F[3] = a * inv;  // (2,2)
    }

    // ------------------------------------------------------------
    // blk_inv_action
    //   L  <- sum_{j=0}^{p-1} alpha[j] * F_j * Es * F_j
    // where F_j = (I + beta[j]*(S-I))^{-1}, and
    //   S = diag(R(theta[0]), ..., R(theta[m-1]), 1?)
    // n = 2*m  or  n = 2*m+1
    // Es, L are n-by-n in column-major storage.
    // ------------------------------------------------------------
    void blk_inv_action(REAL_TYPE *L,
                        const REAL_TYPE *Es,
                        const REAL_TYPE *theta,
                        const REAL_TYPE *alpha,
                        const REAL_TYPE *beta)
    {
        INTE_TYPE n = d;
        INTE_TYPE nblocks = d / 2; // number of 2x2 rotation blocks
        const bool has_scalar = (n == 2 * nblocks + 1);
        const INTE_TYPE sidx = has_scalar ? (2 * nblocks) : -1;

        // Zero out L
        for (INTE_TYPE i = 0; i < n * n; ++i)
            L[i] = 0.0;

        // Workspace: T will hold T = F_j * Es (n x n)
        Work.Zero();
        REAL_TYPE *T = Work.v; // Work is MATX_TYPE(d, d+1), enough for n*n

        // Precompute all 2x2 blocks of F_j,i:
        // F_{j,i} = (I + beta_j (S_i - I))^{-1}, where S_i = R(theta[i])
        // Stored as column-major 2x2 blocks in Fblocks[(j*nblocks + i)*4 .. +3]
        REAL_TYPE *Fblocks = new REAL_TYPE[m * nblocks * 4];

        for (INTE_TYPE j = 0; j < m; ++j)
        {
            REAL_TYPE bet = beta[j];
            for (INTE_TYPE i = 0; i < nblocks; ++i)
            {
                REAL_TYPE phi = theta[i]; // scaled angle for block i
                REAL_TYPE *Fji = &Fblocks[(j * nblocks + i) * 4];
                compute_F_block(Fji, phi, bet);
            }
        }

        // printMat("F blocks (vec4)", Fblocks, 4, m * nblocks);

        // Temporary buffer for right-multiplication (n x 2 panel)
        REAL_TYPE *Ycol = new REAL_TYPE[n * 2];

        // For each Padé node j:
        for (INTE_TYPE j = 0; j < m; ++j)
        {
            const REAL_TYPE aj = alpha[j];

            // ------------------------------------------------------
            // 1) Left action: T = F_j * Es
            //    F_j is block diagonal: each 2x2 block multiplies
            //    the corresponding 2×n row panel of Es.
            // ------------------------------------------------------
            // First, handle all 2x2 rotation blocks (rows)
            for (INTE_TYPE p = 0; p < nblocks; ++p)
            {
                INTE_TYPE r0 = 2 * p;
                INTE_TYPE r1 = r0 + 1;

                const REAL_TYPE *Fjp = &Fblocks[(j * nblocks + p) * 4]; // 2x2
                const REAL_TYPE *Es_panel = &Es[r0];                    // rows r0:r1 of Es, all cols
                REAL_TYPE *T_panel = &T[r0];                            // rows r0:r1 of T, all cols

                // T_panel (2 x n) = Fjp (2 x 2) * Es_panel (2 x n, lda = n)
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            2, n, 2,
                            1.0, Fjp, 2,
                            Es_panel, n,
                            0.0, T_panel, n);
            }

            // printMat("Left action (I + beta (K - I))^{-1} E_s", T, d, d);

            // Scalar row (if present) is unchanged by F_j
            if (has_scalar)
            {
                for (INTE_TYPE col = 0; col < n; ++col)
                    T[sidx + col * n] = Es[sidx + col * n];
            }

            // ------------------------------------------------------
            // 2) Right action: Y_j = T * F_j
            //    Each 2x2 block multiplies the corresponding n×2
            //    column panel of T. We accumulate aj * Y_j into L.
            // ------------------------------------------------------
            // Rotation–rotation column blocks
            for (INTE_TYPE q = 0; q < nblocks; ++q)
            {
                INTE_TYPE c0 = 2 * q;
                INTE_TYPE c1 = c0 + 1;

                const REAL_TYPE *Fjq = &Fblocks[(j * nblocks + q) * 4]; // 2x2
                const REAL_TYPE *T_colpanel = &T[c0 * n];               // n x 2
                REAL_TYPE *Y_panel = Ycol;                              // n x 2

                // Y_panel (n x 2) = T_colpanel (n x 2) * Fjq (2 x 2)
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            n, 2, 2,
                            1.0, T_colpanel, n,
                            Fjq, 2,
                            0.0, Y_panel, n);

                // Accumulate into L: L[:, c0:c1] += aj * Y_panel
                for (INTE_TYPE row = 0; row < n; ++row)
                {
                    L[row + c0 * n] += aj * Y_panel[row + 0 * n];
                    L[row + c1 * n] += aj * Y_panel[row + 1 * n];
                }
            }

            // Scalar column (if present): F_j acts as 1 on it,
            // so Y_j(:, sidx) = T(:, sidx)
            if (has_scalar)
            {
                for (INTE_TYPE row = 0; row < n; ++row)
                    L[row + sidx * n] += aj * T[row + sidx * n];
            }
        }

        delete[] Ycol;
        delete[] Fblocks;
    }

    void _set_R(REAL_TYPE *MatQ, INTE_TYPE ldq) { R.Assign(MatQ, ldq); };
    void _set_A(REAL_TYPE *VecA) { A.Assign(VecA); };

    void Parameter(REAL_TYPE *MatQ, INTE_TYPE ldq, REAL_TYPE *VecA)
    {
        _set_R(MatQ, ldq);
        _set_A(VecA);
    }
    void Parameter(const SkewSchurFactor &SSF) { this->Parameter(SSF.R.v, d, SSF.A.v); };
    void Parameter(ColMat<REAL_TYPE> &MatQ, ArrVec<REAL_TYPE> &VecA) { this->Parameter(MatQ.v, d, VecA.v); };

    void Dlog(REAL_TYPE *MatE, INTE_TYPE ldm)
    {
        // MatE is the perturbation E in the original coordinates.
        // Goal: Y.v will hold dlog(M)[E] in original coordinates at the end.

        // ---------------------------------------------------------
        // 0. Move E into X and transform to Schur coordinates:
        //    X := R^T * E * R
        // ---------------------------------------------------------
        Assign(MatE, ldm); // X := E  (column-major d×d)

        // printMat("R", R.v, d, d);

        // printVec("A", A.v, d / 2);

        // printMat("E", X.v, d, d);

        // Work := R^T * X
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    d, d, d,
                    1.0, R.v, d, X.v, d,
                    0.0, Work.v, d);

        // X := Work * R = R^T * E * R
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    d, d, d,
                    1.0, Work.v, d, R.v, d,
                    0.0, X.v, d);

        // Now X holds E_0 in the Schur frame.

        // printMat("R'ER", X.v, d, d);

        // ---------------------------------------------------------
        // 1. Lyapunov chain for scaling:
        //    For k = 1..s, solve A_k X_k + X_k A_k = X_{k-1}
        //    where A_k = exp(D / 2^k) has angles phi_i = theta_i / 2^k.
        // ---------------------------------------------------------
        const INTE_TYPE nblocks = d / 2; // floor(d/2); for d odd, A.size() = d/2
        VECT_TYPE theta_scaled(nblocks); // local scaled angles phi_i

        if (s > 0)
        {
            for (INTE_TYPE k = 1; k <= s; ++k)
            {
                // phi_i = theta_i / 2^k
                REAL_TYPE inv_scale = std::pow((REAL_TYPE)2.0, (REAL_TYPE)(-k));
                for (INTE_TYPE i = 0; i < nblocks; ++i)
                    theta_scaled[i] = A[i] * inv_scale;

                // Y := X_k solves A_k Y + Y A_k = X  (X = X_{k-1}, Y = X_k)
                blk_lyap(Y.v, X.v, theta_scaled.v);

                // Move result back to X for the next iteration
                // X <- Y
                std::memcpy(X.v, Y.v, sizeof(REAL_TYPE) * d * d);
            }
        }
        // After the loop (or if s==0), X holds E_s at S = exp(D / 2^s).

        // printMat("E_s", X.v, d, d);

        // ---------------------------------------------------------
        // 2. Padé inverse action at S:
        //    L = sum_j alpha[j] (I + beta[j]K)^{-1} E_s (I + beta[j]K)^{-1}
        //    where K = S - I and S has angles phi_i = theta_i / 2^s.
        // ---------------------------------------------------------
        {
            REAL_TYPE inv_scale_s = (s > 0)
                                        ? std::pow((REAL_TYPE)2.0, (REAL_TYPE)(-s))
                                        : (REAL_TYPE)1.0;

            for (INTE_TYPE i = 0; i < nblocks; ++i)
                theta_scaled[i] = A[i] * inv_scale_s;

            // Y := L_log(S)[E_s] via block Padé inverse action
            blk_inv_action(Y.v, // output L
                           X.v, // input E_s
                           theta_scaled.v,
                           alpha, beta);
        }

        // printMat("L", Y.v, d, d);

        // ---------------------------------------------------------
        // 3. Scale back: log(A) = 2^s log(S)  =>  dlog(A)[E] = 2^s L
        // ---------------------------------------------------------
        if (s > 0)
        {
            REAL_TYPE two_to_s = std::pow((REAL_TYPE)2.0, (REAL_TYPE)s);
            cblas_dscal(d * d, two_to_s, Y.v, 1);
        }

        // printMat("2^s L", Y.v, d, d);

        // ---------------------------------------------------------
        // 4. Undo Schur transform: Y := R * Y * R^T
        // ---------------------------------------------------------
        // Work := R * Y
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    d, d, d,
                    1.0, R.v, d, Y.v, d,
                    0.0, Work.v, d);

        // Y := Work * R^T
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    d, d, d,
                    1.0, Work.v, d, R.v, d,
                    0.0, Y.v, d);

        // printMat("2^s RLR'", Y.v, d, d);
        // Now Y.v contains dlog(M)[E] in the original coordinates.
    }

    void Dlog(REAL_TYPE *MatN, INTE_TYPE ldn, REAL_TYPE *MatM, INTE_TYPE ldm)
    {
        Dlog(MatM, ldm);
        Y.Copyto(MatN, ldn);
    }
};