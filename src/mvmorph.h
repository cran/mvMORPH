/* mvmorph.h 2015-01-01 */
/* Julien Clavel        */
#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <Rconfig.h>
#include <Rdefines.h>
#include <R_ext/RS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#define down(x,y) ((x) & ~((y)-1))
#define square(x) (data[(x)]*data[(x)])

La_extern void F77_NAME(dtrttf)(const char* transr, const char* uplo, const La_INT *n,
         const double* a, const La_INT *lda,
         double* arf, La_INT *info FCLEN FCLEN);

La_extern void F77_NAME(dpftrf)(const char* transr, const char* uplo, const La_INT* n,
         double* a, La_INT* info FCLEN FCLEN);

La_extern void F77_NAME(dtfsm)(const char* transr, const char* side, const char* uplo, const char* trans, const char* diag,
        const La_INT* m, const La_INT* n, const double* alpha, const double* a,
        double* b, const La_INT* ldb FCLEN FCLEN FCLEN FCLEN FCLEN);


