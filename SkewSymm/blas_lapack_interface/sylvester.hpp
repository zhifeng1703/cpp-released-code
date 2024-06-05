#pragma once

#include "mkl.h"
#include "type_convention.hpp"
#include "colMajMat.hpp"
#include "schCanFac.hpp"

#include "dgemm.hpp"
#include "dtrsyl.hpp"

INTE_TYPE symsyl(INTE_TYPE sgn, ColMat<REAL_TYPE> &X, const ColMat<REAL_TYPE> &A);
INTE_TYPE symsyl(INTE_TYPE sgn, ColMat<REAL_TYPE> &X, const ColMat<REAL_TYPE> &A, SchurCanonicalFactor &SCF, ColMat<REAL_TYPE> &W);
INTE_TYPE symsyl(INTE_TYPE sgn, ColMat<REAL_TYPE> &X, const ColMat<REAL_TYPE> &A, const ColMat<REAL_TYPE> &C, SchurCanonicalFactor &SCF, ColMat<REAL_TYPE> &W);
INTE_TYPE symsyl(INTE_TYPE sgn, const View_ColMat<REAL_TYPE> &X, const View_ColMat<REAL_TYPE> &A, const View_ColMat<REAL_TYPE> &C, SchurCanonicalFactor &SCF, const View_ColMat<REAL_TYPE> &W);