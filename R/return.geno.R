#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: return.geno.R                                                 #
# Contains: return.geno                                               #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2009, Gabriel R A Margarido                           #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# Function to create diplotypes based on segregation type and linkage phase
return.geno <-
function(segr.type, link.phases) {
  switch(EXPR=segr.type,
         'A.1' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","c","d")),
                  '1.-1'  = return(c("a","b","d","c")),
                  '-1.1'  = return(c("b","a","c","d")),
                  '-1.-1' = return(c("b","a","d","c"))
                  )
         },
         'A.2' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","a","c")),
                  '1.-1'  = return(c("a","b","c","a")),
                  '-1.1'  = return(c("b","a","a","c")),
                  '-1.-1' = return(c("b","a","c","a"))
                  )
         },
         'A.3' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","c","o")),
                  '1.-1'  = return(c("a","b","o","c")),
                  '-1.1'  = return(c("b","a","c","o")),
                  '-1.-1' = return(c("b","a","o","c"))
                  )
         },
         'A.4' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","o","b","o")),
                  '1.-1'  = return(c("a","o","o","b")),
                  '-1.1'  = return(c("o","a","b","o")),
                  '-1.-1' = return(c("o","a","o","b"))
                  )
         },
         'B1.5' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","a","o")),
                  '1.-1'  = return(c("a","b","o","a")),
                  '-1.1'  = return(c("b","a","a","o")),
                  '-1.-1' = return(c("b","a","o","a"))
                  )
         },
         'B2.6' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","o","a","b")),
                  '1.-1'  = return(c("a","o","b","a")),
                  '-1.1'  = return(c("o","a","a","b")),
                  '-1.-1' = return(c("o","a","b","a"))
                  )
         },
         'B3.7' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","a","b")),
                  '1.-1'  = return(c("a","b","b","a")),
                  '-1.1'  = return(c("b","a","a","b")),
                  '-1.-1' = return(c("b","a","b","a"))
                  )
         },
         'C.8' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","o","a","o")),
                  '1.-1'  = return(c("a","o","o","a")),
                  '-1.1'  = return(c("o","a","a","o")),
                  '-1.-1' = return(c("o","a","o","a"))
                  )
         },
         'D1.9' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","c","c")),
                  '1.-1'  = return(c("a","b","c","c")),
                  '-1.1'  = return(c("b","a","c","c")),
                  '-1.-1' = return(c("b","a","c","c"))
                  )
         },
         'D1.10' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","a","a")),
                  '1.-1'  = return(c("a","b","a","a")),
                  '-1.1'  = return(c("b","a","a","a")),
                  '-1.-1' = return(c("b","a","a","a"))
                  )
         },
         'D1.11' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","b","o","o")),
                  '1.-1'  = return(c("a","b","o","o")),
                  '-1.1'  = return(c("b","a","o","o")),
                  '-1.-1' = return(c("b","a","o","o"))
                  )
         },
         'D1.12' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("b","o","a","a")),
                  '1.-1'  = return(c("b","o","a","a")),
                  '-1.1'  = return(c("o","b","a","a")),
                  '-1.-1' = return(c("o","b","a","a"))
                  )
         },
         'D1.13' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","o","o","o")),
                  '1.-1'  = return(c("a","o","o","o")),
                  '-1.1'  = return(c("o","a","o","o")),
                  '-1.-1' = return(c("o","a","o","o"))
                  )
         },
         'D2.14' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("c","c","a","b")),
                  '1.-1'  = return(c("c","c","b","a")),
                  '-1.1'  = return(c("c","c","a","b")),
                  '-1.-1' = return(c("c","c","b","a"))
                  )
         },
         'D2.15' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","a","a","b")),
                  '1.-1'  = return(c("a","a","b","a")),
                  '-1.1'  = return(c("a","a","a","b")),
                  '-1.-1' = return(c("a","a","b","a"))
                  )
         },
         'D2.16' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("o","o","a","b")),
                  '1.-1'  = return(c("o","o","b","a")),
                  '-1.1'  = return(c("o","o","a","b")),
                  '-1.-1' = return(c("o","o","b","a"))
                  )
         },
         'D2.17' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("a","a","b","o")),
                  '1.-1'  = return(c("a","a","o","b")),
                  '-1.1'  = return(c("a","a","b","o")),
                  '-1.-1' = return(c("a","a","o","b"))
                  )
         },
         'D2.18' = {
           switch(EXPR=link.phases,
                  '1.1'   = return(c("o","o","a","o")),
                  '1.-1'  = return(c("o","o","o","a")),
                  '-1.1'  = return(c("o","o","a","o")),
                  '-1.-1' = return(c("o","o","o","a"))
                  )
         }
         )
}

# end of file