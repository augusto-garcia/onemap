#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: group_seq.R                                                   ##
## Contains: group_seq, print.group_seq                                ##
##                                                                     ##
## Written by Cristiane Taniguti                                       ##
##                                                                     ##
## First version: 09/08/2017                                           ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################

## Function to assign markers to preexisting linkage groups

##' Assign markers to preexisting linkage groups
##'
##'
##' Identifies linkage groups of markers combining input \code{sequences} objects with
##' unlinked markers from \code{rf_2pts} object. The results from two-point
##' (pairwise) analysis and the \emph{transitive} property of linkage are used for
##' grouping, as \code{group} function.
##'
##' If the arguments specifying thresholds used to group markers, i.e., minimum
##' LOD Score and maximum recombination fraction, are \code{NULL} (default),
##' the values used are those contained in object \code{input.2pts}. If not
##' using \code{NULL}, the new values override the ones in object
##' \code{input.2pts}.

##' @param input.2pts an object of class \code{rf_2pts}.
##' @param seqs a list of objects of class \code{sequence} or the string
##' "CHROM" if there is \code{CHROM} information available in the input
##' data file.
##' @param unlink.mks a object of class \code{sequence} with the number of
##' the markers to be grouped with the preexisting sequences defined by \code{seqs}
##' parameter. Using the string "all", all remaining markers of
##' the \code{rf_2pts} object will be tested.
##' @param repeated logical. If \code{TRUE}, markers grouped in more than
##' one of the sequences are kept in the output sequences. If \code{FALSE},
##' they are removed of the output sequences.
##' @param LOD a (positive) real number used as minimum LOD score
##' (threshold) to declare linkage.
##' @param max.rf a real number (usually smaller than 0.5) used as
##' maximum recombination fraction to declare linkage.
##' @param min_mks integer defining the minimum number of markers that a provided 
##' sequence (seqs or CHROM) should have to be considered a group. 
##' 
##' @return Returns an object of class \code{group_seq}, which is a list
##'     containing the following components: \item{data.name}{name of
##'     the object of class \code{onemap} that contains the raw
##'     data.} \item{twopt}{name of the object of class \code{rf.2ts}
##'     used as input, i.e., containing information used to assign
##'     markers to linkage groups.} \item{mk.names}{marker names,
##'     according to the input file.} \item{input.seqs}{list with the numbers
##'     of the markers in each inputted sequence}  \item{input.unlink.mks}{numbers of
##'     the unlinked markers in inputted sequence} \item{out.seqs}{list with the
##'     numbers of the markers in each outputted sequence} \item{n.unlinked}{number
##'     of markers that remained unlinked} \item{n.repeated}{number of markers which
##'     repeated in more than one group} \item{n.mar}{total number of markers evaluated}
##'     \item{LOD}{minimum LOD Score to declare linkage.} \item{max.rf}{maximum
##'     recombination fraction to declare linkage.} \item{sequences}{list of outputted
##'     sequences} \item{repeated}{list with the number of the markers that are repeated
##'     in each outputted sequence} \item{unlinked}{number of the markers which remained
##'     unlinked}
##'     
##' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
##' @seealso \code{\link[onemap]{make_seq}} and \code{\link[onemap]{group}}
##'
##' @examples
##' \donttest{
##' data(onemap_example_out) # load OneMap's fake dataset for a outcrossing population
##' data(vcf_example_out) # load OneMap's fake dataset from a VCF file for a outcrossing population
##' comb_example <- combine_onemap(onemap_example_out, vcf_example_out) # Combine datasets
##' twopts <- rf_2pts(comb_example)
##'
##' out_CHROM <- group_seq(twopts, seqs="CHROM", repeated=FALSE)
##' out_CHROM
##'
##' seq1 <- make_seq(twopts, c(1,2,3,4,5,25,26))
##' seq2 <- make_seq(twopts, c(8,18))
##' seq3 <- make_seq(twopts, c(4,16,20,21,24,29))
##'
##' out_seqs <- group_seq(twopts, seqs=list(seq1,seq2,seq3))
##' out_seqs
##' }
##'@export
group_seq <- function(input.2pts, seqs= "CHROM", 
                      unlink.mks="all", 
                      repeated = FALSE, 
                      LOD=NULL, 
                      max.rf=NULL, 
                      min_mks = NULL){
  
  ## checking for correct object
  if(!inherits(input.2pts,"rf_2pts")) stop(deparse(substitute(input.2pts)),
                                     " is not an object of class 'rf_2pts'")
  
  if(!inherits(seqs, "list")){
    if(seqs == "CHROM"){
      ## making CHROM sequences
      CHROM <- unique(input.2pts$CHROM)
      CHROM <- CHROM[!is.na(CHROM)]
      names_seqs <- paste0("CHR",CHROM)
      seqs.int <- list()
      for(i in 1:length(CHROM)) seqs.int[[i]] <- make_seq(input.2pts, CHROM[i])
    } else {
      stop("This option is not available in argument seqs")
    }
  } else{
    ## checking for correct object for seqs argument
    seqs.int <- seqs
    if(!inherits(seqs.int,"list")) stop(deparse(substitute(seqs)),
                                  " is not an object of class 'list'")
    trueseqs <- vector()
    for(i in 1:length(seqs.int)) trueseqs[i] <- inherits(seqs.int[[i]],"sequence")
    if(!all(trueseqs)) stop(" the objects inside the list ",
                            deparse(substitute(seqs)), " are not of class 'sequence'")
    if(is.null(names(seqs.int))) {names_seqs <- paste0("seq",1:length(seqs.int))
    } else { names_seqs <- names(seqs.int)}
  }
  
  ## determining thresholds
  if (is.null(LOD))
    LOD <- input.2pts$LOD
  if (is.null(max.rf))
    max.rf <- input.2pts$max.rf
  
  ## Defining the makers to be tested
  mk_seqs <- sapply(seqs.int, '[[',1)
  
  if(!is.null(min_mks)){
    rm_group <- which(sapply(mk_seqs, length) < min_mks)
    seqs.int[rm_group] <- NULL
    names_seqs <- names_seqs[-rm_group]
  }
  
  mk_seqs <- sapply(seqs.int, '[[',1)
  mk_seqs <- do.call(c, mk_seqs)
  mk_rest <- c(1:input.2pts$n.mar)[-mk_seqs]
  
  if(unlink.mks[1] == "all"){
  } else {
    ## checking for correct object for unlinked.mks argument
    if (!inherits(unlink.mks,"sequence")) {
      stop(" the objects", deparse(substitute(unlink.mks)), " are not of class 'sequence'")
    } else {
      mk_rest <- mk_rest[match(unlink.mks$seq.num, mk_rest)]
      mk_rest <- mk_rest[!is.na(mk_rest)]
    }
  }
  
  if(length(mk_rest) == 0) stop("All markers already have chromosome information. Check min_mks argument.")
  
  ## Grouping
  groups <- new_seqs <- select_group <- seqs_groups <- list()
  same <- vector()
  for(i in 1:length(seqs.int)){
    groups[[i]] <- group(make_seq(input.2pts,c(seqs.int[[i]]$seq.num,mk_rest)),
                         LOD = LOD ,max.rf = max.rf)
    seqs_groups[[i]] <- groups[[i]]$groups[1:length(seqs.int[[i]]$seq.num)]
    same[i] <- length(unique(seqs_groups[[i]])) == 1
    select_group[[i]] <- as.numeric(names(which.max(table(seqs_groups[[i]]))))
    new_seqs[[i]] <- make_seq(groups[[i]], select_group[[i]])
  }
  
  
  if(!all(same)) message(cat("One or more of the provided marker sequences from",deparse(substitute(seqs)),
                     "do not form single linkage groups. The group with the highest number of markers belonging to the sequence will be considered."))
  
  # Find repeated markers
  mks_new_seqs <- lapply(new_seqs, '[[',1)
  all_grouped_mk <- do.call(c, mks_new_seqs)
  repeated_mks <- unique(all_grouped_mk[duplicated(all_grouped_mk)])
  pos_repeated <- lapply(mks_new_seqs,function(x) which(x %in% repeated_mks))
  
  # Unlinked markers
  all <- c(mk_seqs,mk_rest)
  unlinked <- all[-match(all_grouped_mk,all)]
  if(identical(unlinked, integer(0))) unlinked <- NA
  
  names(new_seqs) <- names_seqs
  mk_names <- colnames(input.2pts$data.name$geno)
  
  if(!(identical(repeated_mks, integer(0)) || identical(repeated_mks, numeric(0)))) {
    message("There are one or more markers that grouped in more than one sequence")
    
    # List with repeated markers
    repeated_mks_list <- pos_repeated
    for(i in 1:length(seqs.int)) {
      for(j in 1:length(pos_repeated[[i]]))
        repeated_mks_list[[i]][j] <- mks_new_seqs[[i]][pos_repeated[[i]][j]]
    }
    names(repeated_mks_list) <- names_seqs
    repeated_mks_list[is.na(repeated_mks_list)] <- NULL
    
    # Including or not the repeated in the sequences
    if(repeated){
      structure(list(data.name= input.2pts$data.name, 
                     twopt=input.2pts,
                     mk.names = mk_names, 
                     input.seqs= sapply(seqs.int, '[[',1), 
                     input.unlink.mks= mk_rest,
                     out.seqs = mks_new_seqs, 
                     n.unlinked = length(unlinked[!is.na(unlinked)]),
                     n.repeated = length(unique(unlist(repeated_mks_list))), 
                     n.mar=length(all), 
                     LOD=LOD, 
                     max.rf=max.rf,
                     sequences=new_seqs, 
                     repeated=repeated_mks_list,
                     unlinked= unlinked), class = "group_seq")
    } else {
      new_seqs_unique_temp <- new_seqs_unique <- list()
      for(i in 1:length(seqs.int)) {
        if(identical(pos_repeated[[i]], integer(0))) {
          new_seqs_unique[[i]] <- new_seqs[[i]]
        } else {
          new_seqs_unique_temp[[i]] <- mks_new_seqs[[i]][-pos_repeated[[i]]]
          new_seqs_unique[[i]] <- make_seq(input.2pts, new_seqs_unique_temp[[i]])
          new_seqs_unique[[i]]$twopt <- input.2pts}
      }
      names(new_seqs_unique) <- names_seqs
      structure(list(data.name= input.2pts$data.name, 
                     twopt=input.2pts,
                     mk.names = mk_names, 
                     input.seqs= sapply(seqs.int, '[[',1), 
                     input.unlink.mks= mk_rest,
                     out.seqs = sapply(new_seqs_unique, '[[',1), 
                     n.unlinked = length(unlinked[!is.na(unlinked)]),
                     n.repeated = length(unique(unlist(repeated_mks_list))), 
                     n.mar=length(all), 
                     OD=LOD, 
                     ax.rf=max.rf,
                     sequences=new_seqs_unique, 
                     repeated=repeated_mks_list,
                     unlinked= unlinked), class = "group_seq")}
    
  } else {
    structure(list(data.name= input.2pts$data.name, 
                   twopt=input.2pts,
                   mk.names = mk_names, 
                   input.seqs= sapply(seqs.int, '[[',1), 
                   input.unlink.mks= mk_rest,
                   out.seqs = mks_new_seqs, 
                   n.unlinked = length(unlinked[!is.na(unlinked)]),
                   n.repeated = 0, 
                   n.mar=length(all), 
                   LOD=LOD, 
                   max.rf=max.rf,
                   sequences=new_seqs, 
                   repeated=NA,
                   unlinked= unlinked), class = "group_seq")
  }
}

##' Show the results of grouping markers to preexisting sequence
##'
##' It shows the groups sequences, the repeated markers, as well as the unlinked markers.
##'
##' @aliases print.group_seq
##' @param x an object of class group_seq
##'
##' @param detailed logical. If \code{TRUE} the markers in each
##'     linkage group sequence are printed.
##'
##' @param ... currently ignored
##'
##' @return No return value, called for side effects
##' @keywords internal
##' 
##' @method print group_seq
##' 
##' @export
print.group_seq <- function(x, detailed=TRUE,...) {
  
  ## checking for correct object
  if(!inherits(x,"group_seq")) stop(deparse(substitute(x))," is not an object of class 'group_seq'")
  
  cat("  This is an object of class 'group_seq'\n")
  
  ## criteria
  cat("  Criteria used to assign markers to groups:\n")
  cat("    LOD =", x$LOD, ", Maximum recombination fraction =",
      x$max.rf, "\n")
  
  ## printing summary
  cat("\n  No. markers in input sequences:\n")
  
  for(i in 1:length(x$sequences)){ cat("                      ",names(x$sequences)[[i]],":  ",
                                       length(x$input.seqs[[i]]), "markers\n")}
  cat("\n  No. unlinked input markers:  ", length(x$input.unlink.mks), "markers\n")
  cat("\n  No. markers in output sequences:\n")
  for(i in 1:length(x$sequences)){ cat("                      ",names(x$sequences)[[i]],":  ",
                                       length(x$out.seqs[[i]]), "markers\n")}
  cat("  No. unlinked:                ", x$n.unlinked, "markers\n")
  cat("  No. repeated:                ", x$n.repeated, "markers\n")
  
  if (detailed) {
    ## printing detailed results (markers in each linkage group)
    cat("\n  Printing output sequences:")
    for (i in 1:length(x$sequences)) {
      cat("\n  Group", names(x$sequences)[[i]], ":", length(x$out.seqs[[i]]) , "markers\n    ")
      cat(x$mk.names[x$sequences[[i]]$seq.num], "\n")
    }
    cat("\n  Unlinked markers:", x$n.unlinked ,"markers\n    ")
    if(x$n.unlinked==0) {cat("\n")} else cat(x$mk.names[x$unlinked], "\n")
    cat("\n  Repeated markers:", x$n.repeated, " markers\n     ")
    if(x$n.repeated==0) {
      cat("\n")
    } else {
      for (i in 1:length(x$repeated)) {
        cat("\n  Group", names(x$repeated)[[i]], ":", length(x$repeated[[i]]) , "markers\n    ")
        cat(x$mk.names[x$repeated[[i]]], "\n")
      }
    }
  }
}
