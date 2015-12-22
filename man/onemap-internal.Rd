\name{onemap-internal}
\alias{acum}
\alias{adjust.rf.ril}
\alias{avoid.reverse}
\alias{codif.data}
\alias{comb}
\alias{comb.ger}
\alias{cr2pts}
\alias{cr3pts}
\alias{Cindex}
\alias{diplo}
\alias{draw.order}
\alias{draw.rf.inter}
\alias{draw.try}
\alias{est.map.c}
\alias{est_map_hmm_f2}
\alias{est_rf_bc}
\alias{est_rf_f2}
\alias{est_rf_out}
\alias{get.bins}
\alias{get_mat_rf_in}
\alias{get_vec_rf_in}
\alias{perm.pars}
\alias{perm.tot}
\alias{phases}
\alias{plotFunction.out}
\alias{rem.amb.ph}
\alias{return.geno}
\alias{ser.ord}
\alias{print.bc.onemap}
\alias{print.compare}
\alias{print.f2.onemap}
\alias{print.order}
\alias{print.outcross}
\alias{print.rf.3pts}
\alias{print.riself.onemap}
\alias{print.risib.onemap}
\alias{print.sequence}
\title{Internal onemap functions}
\description{
  Functions from \pkg{onemap} not to be directly called by the user.
}
\usage{
acum(w)
avoid.reverse(input)
adjust.rf.ril(r,type,expand)
codif.data(geno.in,segr.type.in)
comb(x,y)
comb.ger(f)
cr2pts(mrk1, mrk2, segr.type1, segr.type2)
cr3pts(mrk1, mrk2, mrk3, segr.type1, segr.type2, segr.type3)
Cindex(order,r)
diplo(w, seq.num, seq.phases)
draw.order(map.input)
draw.rf.inter(info.graph, scale, cex)
draw.try(base.input, try.input, pos=NULL)
est.map.c(geno,type,phase,rec,verbose,tol=1e-2)
est_rf_out(geno, mrk, seg_type, nind, verbose)
est_rf_f2(geno, mrk, seg_type, nind, verbose)
est_rf_bc(geno, mrk,  nind, type, verbose)
est_map_hmm_f2(geno, rf.vec, verbose)
get.bins(geno, exact)
get_mat_rf_in(input.seq)
get_vec_rf_in(input.seq)
perm.pars(v)
perm.tot(v)
phases(input.seq, LOD=0, max.rf=0.50)
plotFunction.out(info.graph, cex, main)
rem.amb.ph(M,w,seq.num)
return.geno(segr.type, link.phases)
ser.ord(r,i)
\method{print}{bc.onemap}(x, \dots)
\method{print}{compare}(x, \dots)
\method{print}{f2.onemap}(x, \dots)
\method{print}{order}(x, \dots)
\method{print}{outcross}(x, \dots)
\method{print}{rf.3pts}(x, \dots)
\method{print}{sequence}(x, \dots)
\method{print}{riself.onemap} (x, \dots)
\method{print}{risib.onemap} (x, \dots)
}
\author{Gabriel Rodrigues Alves Margarido and Marcelo Mollinari}
\keyword{internal}
