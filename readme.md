# Instruction of the `dexp` Project

Source Code and Experiments for
“The Exponential of Skew-Symmetric Matrices: A Nearby Inverse and Efficient Computation of Derivatives”
Authors: Zhifeng Deng, P.-A. Absil, Kyle Gallivan, Wen Huang
arXiv: https://arxiv.org/abs/2506.18302

This repository contains the C++ source code, MATLAB interfaces, and experiment drivers used to generate the numerical results in Figures 3–8 of the paper.

---

## Folder Structure

```
dexp/
│
├── main.cpp          # C++ executable driver
├── basics.hpp        # Basic experiment utilities
│
├── skewblas/         # C++ source code: Schur blocks, expm, dexp, DK formulas, etc.
│
├── test/             # C++ experiment codes (Figures 3–8)
│
├── cpp-api/          # C++ dynamic library settings
│
├── matlab-mex/       # MATLAB MEX gateway to the C++ utilities
│
├── matlab-dll/       # MATLAB callable C++ dynamic library (Windows only)
│
└── figures/          # Figures inputs and outputs of the experiments.


```

---

## 1. Run codes in MATLAB

These codes are ran and tested under the MATLAB 2024a under the following environment:

- CPU/RAM: 13th Gen Intel(R) Core(TM) i5-13600KF, 14 Cores, 3.50 GHz / 16 GB. 731
- OS: Windows 11
- BLAS and LAPACK Library: Intel(R) Math Kernel Library.

### 1.1 MATLAB MEX Gateway (dexp/matlab-mex)

To reproduce experiments in Figures 5(a), 7(a), and 8:

       cd dexp/matlab-mex
       test_dexp
       test_expm

The scripts will automatically compile the MEX files when needed.

> **Warning:** Run multiple times for stable timing.

Available MEX utilities:

- Schur decomposition: mex_sblas_ss2schur, mex_sblas_so2schur
- Matrix exponential: mex_sblas_expm, mex_pade_expm
- Differential of expm: sblas_dexp, pade_dexp, dk_dexp

### 1.2 MATLAB + Windows DLL (dexp/matlab-dll)

The MATLAB-callable dll is compiled under the following environment:

- CPU/RAM: 13th Gen Intel(R) Core(TM) i5-13600KF, 14 Cores, 3.50 GHz / 16 GB. 731
- OS: Windows 11
- BLAS and LAPACK Library: Intel(R) Math Kernel Library.

To use the MATLAB-callable DLL:

       cd dexp/matlab-dll
       test_dexp
       test_expm

> **Warning:** MATLAB may crash when repeatedly passing large matrices
> (100×100 or larger) between MATLAB and the DLL due to unstable
> MATLAB memory-management behavior.

---

## 2. Run codes in C++

These codes are ran and tested under the following environment:

- CPU/RAM: 13th Gen Intel(R) Core(TM) i5-13600KF, 14 Cores, 3.50 GHz / 16 GB. 731
- OS: Windows 11, Windows Subsystem Linux (Ubuntu 22.04.5 LTS). 732
- BLAS and LAPACK Library: Intel(R) Math Kernel Library.
- Requirements: - g++ (C++17 or later) - Intel OneMKL (BLAS/LAPACK):
  https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html

The compilation relies on the environment variable MKLROOT, which specifies the location of Intel OneAPI MKL (typically /path/to/your/intel/oneapi/mkl/latest). This variable is automatically set when the Intel oneAPI environment is initialized, e.g.
```bash
source /path/to/your/intel/oneapi/setvars.sh
```
If MKLROOT is not defined, the Makefile contains the fallback line at line 14
```bash
MKLROOT ?= /opt/intel/oneapi/mkl/latest
```
which uses the default installation path. Users may edit this line if MKL is installed elsewhere, but this is not necessary when setvars.sh has been sourced correctly.

---

### 2.1 Reproducing Paper Experiments

For reproducing paper experiments, compile the test codes /test/*_test.cpp with the line:
```bash
make tests
```
which produces two executables: elapsed_test and error_test. For Figure 3(a), 4-7, 8(a) in the paper: run
```bash
./elapsed_test -task forward
```
for the differentiation, and run ./elapsed_test -task inverse 
```bash
./elapsed_test -task inverse
```
for the inversed differentiation. For Figures 3(b), 8(b) in the paper: run     ```bash
./error_test
```

---
