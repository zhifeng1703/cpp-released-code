// #include "mkl.fi"
#include "type_convention.hpp"

// extern "C" void dlaqtr(BOOL_TYPE ltran, BOOL_TYPE lreal, INTE_TYPE n, REAL_TYPE *t, INTE_TYPE ldt, REAL_TYPE *b, REAL_TYPE *w, REAL_TYPE scale, REAL_TYPE *x, REAL_TYPE *work, INTE_TYPE *info);

INTE_TYPE my_dlaqtr_real(BOOL_TYPE trans, INTE_TYPE n, REAL_TYPE *T, INTE_TYPE ldt, REAL_TYPE scale, REAL_TYPE *x, REAL_TYPE *work);
