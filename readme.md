`dexp` – Source Code and Experiments for
“The Exponential of Skew-Symmetric Matrices: A Nearby Inverse and Efficient Computation of Derivatives”
Authors: Zhifeng Deng, P.-A. Absil, Kyle Gallivan, Wen Huang
arXiv: https://arxiv.org/abs/2506.18302

This repository contains the C++ source code, MATLAB interfaces, and experiment drivers used to generate the numerical results in Figures 3–8 of the paper.

---

## Folder Structure

dexp/
skewblas/ C++ source code: Schur blocks, expm, dexp, DK formulas, utilities.
test/ C++ experiment drivers (Figures 3–8).
matlab-mex/ MATLAB MEX gateway to the C++ utilities.
matlab-dll/ MATLAB callable C++ dynamic library for Windows.

---

## MATLAB Usage

### 1. MATLAB MEX Gateway (dexp/matlab-mex)

To reproduce experiments in Figures 5(a), 7(a), and 8:

       cd dexp/matlab-mex
       test_dexp
       test_expm

The scripts will automatically compile the MEX files when needed.
Run multiple times for stable timing.

Available MEX utilities:
Schur decomposition: mex_sblas_ss2schur, mex_sblas_so2schur
Matrix exponential: mex_sblas_expm, mex_pade_expm
Differential of expm: sblas_dexp, pade_dexp, dk_dexp

### 2. MATLAB + Windows DLL (dexp/matlab-dll)

To use the MATLAB-callable DLL:

       cd dexp/matlab-dll
       test_dexp
       test_expm

Warning: MATLAB may crash when repeatedly passing large matrices
(100×100 or larger) between MATLAB and the DLL due to unstable
MATLAB memory-management behavior.

---

## Compiling the C++ Source Code

Requirements: - g++ (C++17 or later) - Intel OneMKL (BLAS/LAPACK):
https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html

Before compiling, edit the Makefile and set:

    MKLROOT = /path/to/your/oneMKL/installation

To compile the main code:

    make

---

### Reproducing Paper Experiments

Figure 3:
Replace dexp/main.cpp with dexp/test/elapsed_test.cpp
Run: make

Figures 4–8:
Same procedure:
Replace dexp/main.cpp with dexp/test/error_test.cpp
Run: make

---
