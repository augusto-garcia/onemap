#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP est_hmm_bc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP est_hmm_f2(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP est_hmm_out(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP est_rf_bc_wrap(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP est_rf_f2_wrap(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP est_rf_out_wrap(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_bins(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"est_hmm_bc",      (DL_FUNC) &est_hmm_bc,      5},
  {"est_hmm_f2",      (DL_FUNC) &est_hmm_f2,      5},
  {"est_hmm_out",     (DL_FUNC) &est_hmm_out,     7},
  {"est_rf_bc_wrap",  (DL_FUNC) &est_rf_bc_wrap,  5},
  {"est_rf_f2_wrap",  (DL_FUNC) &est_rf_f2_wrap,  5},
  {"est_rf_out_wrap", (DL_FUNC) &est_rf_out_wrap, 5},
  {"get_bins",        (DL_FUNC) &get_bins,        3},
  {NULL, NULL, 0}
};

void R_init_onemap(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}