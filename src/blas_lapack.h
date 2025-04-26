#ifndef HJ_BLAS_LAPACK_H_
#define HJ_BLAS_LAPACK_H_

#ifdef __cplusplus
namespace f2c {  // prevent typedef pollution
extern "C" {
#endif
#include <f2c.h>
#ifdef __cplusplus
} }
#endif

#undef min
#undef max
#undef abs

#ifdef __cplusplus
namespace clapack {  // for using f2c locally
using namespace f2c;
#endif
#ifdef _MSC_VER
#undef small
#endif
// clapack.h already include extern "C"
#include "clapack.h"
#ifdef __cplusplus
}
#endif


// cblas.h already include extern "C"
// #include <cblas.h> // try to not use cblas

#ifdef __cplusplus
extern "C" {
#endif
#pragma push_macro("ADD_")
#define ADD_
#ifndef F77_GLOBAL  
#  define F77_GLOBAL(f, F) f##_ // depends on the specific lib
#endif  
#include "cblas_f77.h"
#pragma pop_macro("ADD_")
#ifdef __cplusplus
}
#endif

#endif
