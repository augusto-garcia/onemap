\name{onemap-internal}
\alias{acum}
\alias{adjust_rf_ril}
\alias{avoid_reverse}
\alias{check_linkage}
\alias{codif_data}
\alias{comb}
\alias{comb_ger}
\alias{compare_inbred_f2}
\alias{compare_inbred_bc}
\alias{compare_outcross}
\alias{Cindex}
\alias{diplo}
\alias{draw_order}
\alias{draw_rf_inter}
\alias{draw.try}
\alias{est_map_hmm_f2}
\alias{est_map_hmm_bc}
\alias{est_rf_bc}
\alias{est_rf_out}
\alias{get_bins}
\alias{get_mat_rf_in}
\alias{get_vec_rf_in}
\alias{get_mat_rf_out}
\alias{get_twopt_info}
\alias{get_vec_rf_out}
\alias{perm_pars}
\alias{perm_tot}
\alias{phases}
\alias{rem_amb_ph}
\alias{return_geno}
\alias{ripple_seq_inbred}
\alias{ripple_seq_outcross}
\alias{ser_ord}
\alias{try_seq_inbred_f2}
\alias{try_seq_inbred_bc}
\alias{try_seq_outcross}
\title{Internal onemap functions}
\description{
  Functions from \pkg{onemap} not to be directly called by the user.
}
\usage{
acum(w)
avoid_reverse(input)
adjust_rf_ril(r,type,expand)
check_linkage(i, s, cl, geno, st, max.rf, LOD)
codif_data(geno.in,segr.type.in,cross = c("outcross", "f2", "backcross",
"riself", "risib"))
comb(x,y)
comb_ger(f)
compare_inbred_bc(input.seq, n.best, tol, verbose = FALSE)
compare_inbred_f2(input.seq, n.best, tol, verbose = FALSE)
compare_outcross(input.seq, n.best, tol, verbose)
Cindex(order,r)
diplo(w, seq.num, seq.phases)
est_rf_out(geno, mrk, seg_type, nind, verbose)
est_rf_bc(geno, mrk,  nind, type, verbose)
est_map_hmm_f2(geno, error, rf.vec, verbose, tol)
est_map_hmm_bc(geno, error, rf.vec, verbose, tol)
get_bins(geno, exact)
get_mat_rf_in(input.seq, LOD, max.rf, min.LOD)
get_twopt_info(twopt, small, big)
get_vec_rf_in(input.seq, LOD, acum)
get_mat_rf_out(input.seq, LOD, max.rf, min.LOD)
get_vec_rf_out(input.seq, LOD, max.rf, min.LOD, acum)
perm_pars(v)
perm_tot(v)
phases(input.seq, LOD=0, max.rf=0.50)
rem_amb_ph(M,w,seq.num)
return_geno(segr.type, link.phases)
ripple_seq_inbred(input.seq, ws, ext.w, LOD, tol, verbose)
ripple_seq_outcross(input.seq, ws, LOD, tol, verbose)
ser_ord(r,i)
try_seq_inbred_f2(input.seq, mrk, tol, pos, verbose)
try_seq_outcross(input.seq, mrk, tol, pos, verbose)
\method{print}{compare}(x, \dots)
\method{print}{onemap}(x, \dots)
\method{print}{order}(x, \dots)
\method{print}{order}(x, \dots)
\method{print}{sequence}(x, \dots)
}
\author{Gabriel Rodrigues Alves Margarido and Marcelo Mollinari}
\keyword{internal}
