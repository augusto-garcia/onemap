#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: map_func.R                                                    #
# Contains: haldane, kosambi, set.map.fun                             #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido and Marcelo Mollinari  #
# copyright (c) 2007-9, Gabriel R A Margarido                         #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 11/07/2007                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################
haldane <-
function(rcmb) {
  # rcmb must be a number between 0 and 0.5
  if (is.numeric(rcmb)) {
    if (rcmb >= 0 && rcmb <= 0.5)
      return(100*-(1/2)*log(1-2*rcmb))  # haldane mapping fuction
    else stop("the recombination fraction must be between 0 and 0.5")
  }
  else stop("the argument to 'haldane' function must be ",
            dQuote("numeric"))
}


kosambi <-
function(rcmb) {
  # rcmb must be a number between 0 and 0.5
  if (is.numeric(rcmb)) {
    if (rcmb >= 0 && rcmb <= 0.5)
      return(100/4*log(abs((1+2*rcmb)/(1-2*rcmb))))  # kosambi mapping fuction
    else stop("the recombination fraction must be between 0 and 0.5")
  }
  else stop("the argument to 'kosambi' function must be ",
            dQuote("numeric"))
}

set.map.fun<-
  function(type=c("kosambi", "haldane")){
    type<-match.arg(type,c("kosambi", "haldane"))
    if(type=="haldane") assign(".map.fun",  "haldane", envir = .onemapEnv)
    else assign(".map.fun",  "kosambi", envir = .onemapEnv)
  }

# end of file
