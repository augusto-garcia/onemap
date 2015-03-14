#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: avoid.reverse.R                                               #
# Contains: avoid.reverse                                             #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2009, Marcelo Mollinari                               #
#                                                                     #
# First version: 01/15/2010                                           #
# Last update: 01/15/2010                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################


## This is to avoid reverse orders
avoid.reverse<-function(input) {
  corr1 <- cor.test(input,1:length(input),method="spearman",alternative="two.sided")$estimate
  corr2 <- cor.test(rev(input),1:length(input),method="spearman",alternative="two.sided")$estimate
  if (corr1 > corr2) output<-input
  else if (corr1 < corr2) output<-rev(input)
  else {
    rand <- sample(2,1)
    if (rand == 1) output<-input else output<-rev(input)
  }
  output
}
