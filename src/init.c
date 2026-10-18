#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP C_lsap(SEXP);
extern SEXP C_perm_dcov(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"C_lsap", (DL_FUNC) &C_lsap, 1},
  {"C_perm_dcov", (DL_FUNC) &C_perm_dcov, 3},
  {NULL, NULL, 0}
};

void R_init_HrankdCov(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
