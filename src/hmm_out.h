#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

RcppExport SEXP est_hmm_out(SEXP geno_R, SEXP type_R, SEXP phase_R, SEXP rf_R, SEXP verbose_R, SEXP tol_R);
double emit_out(int obs_gen, int true_gen, double error_prob, int mark_type);
double step_out(int gen1, int gen2, int phase, double rf);
double nrec_out(int gen1, int gen2, int phase);
