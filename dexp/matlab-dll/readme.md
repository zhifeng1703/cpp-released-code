# Matrix Exponential & Differential Utilities (MATLAB + sBLAS)

This repository contains MATLAB front-end wrappers and experiments for testing and benchmarking fast C++/MKL implementations of:

- expm(S)
- Dexp_S(X)   (forward differential)
- (Dexp_S)^{-1}(Y)   (inverse differential)

All core numerical kernels are implemented in skewblas.dll, and MATLAB interfaces are provided in ./src.

Directory layout:

src/        MATLAB wrapper functions (sblas_expm, sblas_dexp, pade_expm, pade_dexp, dk_dexp, ...)
utest/      Unit tests validating correctness of DLL wrappers
test/       Timing and scaling experiments (matrix exponential and differential)
lib/        skewblas.dll and skewblas_api.h


============================================================
1. Calling the basic utilities (./src/*.m)
============================================================

Matrix exponential:
    Q = sblas_expm(S);          % DLL-based exponential
    Q = pade_expm(S, p);        % Pade approximant with order p (3 or 13)
    Q = expm(S);                % MATLAB built-in (reference)

Input:
    S — real n×n skew-symmetric matrix

Output:
    Q — orthogonal matrix expm(S)


Differential of the exponential (forward):
    Y = sblas_dexp(S, X, 'fwd');
    Y = pade_dexp(S, X, 'fwd', p);
    Y = dk_dexp(S, X, 'fwd');

Differential of the exponential (inverse):
    X = sblas_dexp(S, Y, 'inv');
    X = pade_dexp(S, Y, 'inv', p);
    X = dk_dexp(S, Y, 'inv');

Inputs:
    S — skew-symmetric base point
    X or Y — direction matrix, same size as S
    p — Pade order (optional)

Outputs:
    Forward: Y = Dexp_S(X)
    Inverse: X = (Dexp_S)^{-1}(Y)


============================================================
2. Running the unit tests (./utest/utest_xxx.m)
============================================================

Unit tests validate correctness of the wrapper functions.  
Examples:

    utest_expm
    utest_sblas_dexp
    utest_pade_dexp
    utest_dk_dexp

Inputs:
    Most unit tests take no required inputs.
    Some allow optional arguments such as dimension or scale.

Outputs:
    Printed diagnostic information (timings, errors).
    Many tests also generate convergence plots.
    Most do not return variables unless documented.


============================================================
3. Running the experiments (./test_expm.m and ./test_dexp.m)
============================================================

These scripts benchmark timing and scaling of all implementations.

------------------------------------------------------------
Matrix exponential experiment:
------------------------------------------------------------

    elapsed_times = test_expm(3:50);

Input:
    dims — vector of matrix sizes (default 3:50)

Output:
    elapsed_times — (#dims × 4) matrix of averaged timings:
        1. sBLAS expm
        2. Pade p=3
        3. Pade p=13
        4. MATLAB expm

Plots:
    • Elapsed time vs. dimension
    • Frobenius error vs. dimension (log scale)


------------------------------------------------------------
Differential experiment:
------------------------------------------------------------

    [elapsed_forward, elapsed_inverse] = test_dexp(3:50, scale);

Inputs:
    dims  — vector of matrix dimensions
    scale — scaling of S (default 1)

Outputs:
    elapsed_forward — (#dims × 4)
        1. sBLAS Dexp (fwd)
        2. Pade p=3
        3. Pade p=13
        4. DK (forward)

    elapsed_inverse — (#dims × 4)
        1. sBLAS Dexp (inv)
        2. Pade p=1
        3. Pade p=7
        4. DK (inverse)

Plots:
    • Forward Dexp timing vs dimension
    • Inverse Dexp timing vs dimension
    (Both linked on x-axis)


============================================================
DLL Loading Notes
============================================================

All utilities automatically call:

    loadlibrary('skewblas.dll', 'skewblas_api.h')

if the DLL is not already loaded.

Place skewblas.dll and skewblas_api.h in ./lib.


============================================================
End of README
============================================================
