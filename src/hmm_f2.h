/* Written by Marcelo Mollinari
   Adapted from hmm_main.c, hmm_f2.c and util.c (found in the R package qtl)
   copyright (c) 2001-10, Karl W Broman    

   These codes are under the GNU General Public License, version 3
   A copy of the GPL 3 is available at http://www.r-project.org/Licenses/GPL-3
*/

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
double stepf_f2(int gen1, int gen2, double rf);
double nrecf_f2(int gen1, int gen2);
RcppExport SEXP est_hmm_f2(SEXP geno_R, SEXP rf_R, SEXP verbose_R);
