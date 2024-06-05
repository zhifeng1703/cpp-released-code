This project (cpp-release-code/SkewSymm) compiles tasks assoicated to the set of skew symmetric matrices and spectial orthogonal matrices. It runs on Linux environment with Intel CPU and requries the Intel MKL as BLAS library. 


============================================Environment============================================

# An Intel MKL installment is included in the following path.
#	cplusplus-release-code/env_setup/l_onemkl_p_2024.0.0.49673_offline.sh
# It is too big for GitHub Repo, please download it from the website via command (wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/2f3a5785-1c41-4f65-a2f9-ddf9e0db3ea0/l_onemkl_p_2024.1.0.695_offline.sh)

Use the following terminal command to install the library, or refer to https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html for more options.
	sudo sh path-to-installment/_onemkl_p_2024.1.0.695_offline.sh

The default path of the installed Intel MKL is in (/opt/intel/oneapi/mkl/latest). If the library is installed elsewhere, UPDATE the variable (MKLROOT = /opt/intel/oneapi/mkl/latest) in makefile to the corresponding path.

============================================Compilation============================================

The tasks implemented in this project are stored in the several cpp files in the (cpp-release-code/SkewSymm/test) folder. To complile, simply run (make).

To switch to different task, run (cp test/other-task.cpp main.cpp) before (make).

To enable/disable compiler optimization or debug feasure, change the makefile variable (EXTRA = -o3) currently on the maximal optimized compilation (-o3). Users may use (EXTRA = -o1), (EXTRA = -o2), (EXTRA = -o3) or disable EXTRA for different compilation. To enable the debug mode, use (EXTRA = -g).

============================================Execution and Flags============================================

Various tasks generate the executable file under the same name: skewsymm, but accept different flags. All flags have a default values and therefore the execution runs well with no flag provided. To modified the flags, for example, with 10 for -dim and 9527 for -seed, run (./skewsymm -dim 10 -seed 9527)

dexp_real_cmpx_pade_test.cpp 
	It runs sanity check on Dexp_A[M] = exp(A)N with randomly sampled skew symmetric matrix A and M. The real formula, the complex formula and the Pade algorithm are implemented. It accepts the following keys:
		-dim (specifies the size of A and M).
		-seed (sepcifies seed for the random number generator).

dexp_elapsed_time_test.cpp and dlog_elapsed_time_test.cpp
	They time the execution of M |-> N or N |-> M from Dexp_A[M] = exp(A)N with random skew symmetric matrix A and M with various size. It accepts the following keys:
		-dim (specifies the largest dimension d > 4 of A and M. The routine will time the execution on matrices from size 4 x 4, 5 x 5, ..., to d x d).
		-loop (specifies l the number of execution. The same task will be executed l times and the elapsed time will be recorded as either the minimum or the average among them).
		-file (specifies the path to the binary file that save the recorded elapsed time).
		-seed (sepcifies seed for the random number generator).
	Example command:
		./skewsymm -file ~/output -dim 100 -loop 50 -seed 9527

stiefel_geodesic_elapsed_time_test.cpp
	They time the BCH algorithm and the Newton algorithm on solving St_{n,p} endpoint geodesic problem. It accepts the following keys:
		-row (specifies n in St_{n,p}).
		-col (specifies p in St_{n,p}).
		-iter (specifies the maximum iterations allowed in the iterative algorithms).
		-gmres-restart (specifies the restart dimension of the Krylov space generated in GMRES that is required in the Newton algorithm).
		-gmres-maxiter (specifies the maximum iteration allowed in the GMRES that is required in the Newton algorithm).
		-scale-beg (specifies the smallest separation of the endpoints on St_{n,p}, measured by matrix-2-norm)
		-scale-end (specifies the largest separation of the endpoints on St_{n,p}, measured by matrix-2-norm)
		-scale-num (specifies the number of various separation taken within [-scale-beg, -scale-end].)
		-loop (specifies l the number of execution. The same task will be executed l times and the elapsed time will be recorded as either the minimum or the average among them).
		-file (specifies the path to the binary file that save the recorded elapsed time).
		-seed (sepcifies seed for the random number generator).
	Example command:
		./skewsymm -file ~/output -row 40 -col 20 -iter 100 -loop 50 -scale-beg 1.0 -scale-end 3.5 scale-num 100 -seed 9527 -gmres-restart 20 -gmres-maxiter 200





