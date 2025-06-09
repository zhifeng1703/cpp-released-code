To compile the source code, please install gpp and OneMKL for BLAS and LAPACK support (https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) in your linux environment. Make sure to set the correct path of the keyword ``MKLROOT'' in the makefile.

To compile the experiment code of the paper ``The Exponential of Skew-Symmetric Matrices: A Nearby Inverse and Efficient Computation of Derivatives", run make under the directory dexpSkewSymm.

More source code of different experiments are stored under the text sub-directory and you may overwrite the main.cpp with any of them to compile a different experiment setup.
