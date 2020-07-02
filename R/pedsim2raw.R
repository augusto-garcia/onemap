globalVariables(c("read.table"))

#' Converts the output of PedigreeSim to onemap raw file
#' 
#' 
#' @param cross string defining the cross type "outcross" or "f2 intercross"
#' @param genofile pedigreeSim output .dat
#' @param parent1 string defining the first parent ID
#' @param parent2 string defining the seconde parent ID
#' @param f1 string defining the F1 ID, if cross type "f2 intercross
#' @param out.file string defining the name of the output file
#' @param miss.perc double defining the percentage of missing data to be simulated
#' 
#' 
#' @examples 
#' \dontrun{
#'  # Outcrossing population
#'  pedsim2raw(cross="outcross", genofile = system.file("extdata/sim_out_genotypes.dat", 
#'             package = "onemap"), parent1 = "P1", parent2 = "P2", 
#'              out.file = "sim_out.example.raw")
#' 
#'  df <- read_onemap(inputfile = "sim_out.example.raw")
#' }
#' @export
pedsim2raw <- function(cross = c("outcross", "f2 intercross"),
                       genofile = "sim_out_genotypes.dat",
                       parent1 = "P1",
                       parent2 = "P2",
                       f1 = "F1",
                       out.file = "sim_out.raw",
                       miss.perc = 10
){
  
  genofile <- read.table(genofile, header = T, stringsAsFactors = F)
  n.mk <- dim(genofile)[1]
  
  idxP1.1 <- which(colnames(genofile) == paste0(parent1,"_1"))
  idxP1.2 <- which(colnames(genofile) == paste0(parent1,"_2"))
  idxP2.1 <- which(colnames(genofile) == paste0(parent2,"_1"))
  idxP2.2 <- which(colnames(genofile) == paste0(parent2,"_2"))
  
  if(cross == "f2 intercross"){
    n.ind <- (dim(genofile)[2]-1)/2 -3
    idx.F1.1 <- which(colnames(genofile) == paste0(f1,"_1"))
    idx.F1.2 <- which(colnames(genofile) == paste0(f1,"_2"))
    genodat <- genofile[,-c(1,idxP1.1,idxP1.2,idxP2.1,idxP2.2, idx.F1.1, idx.F1.2)]
  } else{
    n.ind <- (dim(genofile)[2]-1)/2 -2
    genodat <- genofile[,-c(1,idxP1.1,idxP1.2,idxP2.1,idxP2.2)]
  }
  
  # Types
  idx.odd <- seq(from=1, to=(n.ind)*2-1, 2)
  idx.even <- seq(from=2, to=(n.ind)*2, 2)
  
  segr.type.in <- rep(NA, n.mk)
  gt_matrix <- matrix(rep(NA,(n.ind)*n.mk), nrow = n.mk, ncol = n.ind)
  if(cross=="outcross"){
    idx <- which(genofile[,idxP1.1] == "a" & genofile[,idxP1.2] == "b" | 
                   genofile[,idxP1.1] == "b" & genofile[,idxP1.2] == "a")
    
    if(length(idx)>0){
      idx2 <- which(genofile[idx,idxP2.1] == "c" & genofile[idx,idxP2.2] == "d" | 
                      genofile[idx,idxP2.1] == "d" & genofile[idx,idxP2.2] == "c")
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "A.1" #ab x cd
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "d" |
                                      genodat[idx[idx2],idx.odd] == "d" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "ad"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "c" |
                                      genodat[idx[idx2],idx.odd] == "c" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "ac"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "c" |
                                      genodat[idx[idx2],idx.odd] == "c" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "bc"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "d" |
                                      genodat[idx[idx2],idx.odd] == "d" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "bd"
      }
      
      idx2 <- which(genofile[idx,idxP2.1] == "a" & genofile[idx,idxP2.2] == "c" |
                      genofile[idx,idxP2.1] == "c" & genofile[idx,idxP2.2] == "a")
      if(length(idx2)>0){
        segr.type.in[idx[idx2]] <- "A.2" # ab x ac
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "a" ,arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "c" |
                                      genodat[idx[idx2],idx.odd] == "c" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "ac"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "a" |
                                      genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "ba"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "c" |
                                      genodat[idx[idx2],idx.odd] == "c" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "bc"
      }
      
      idx2 <- which(genofile[idx,idxP2.1] == "c" & genofile[idx,idxP2.2] == "o" |
                      genofile[idx,idxP2.1] == "o" & genofile[idx,idxP2.2] == "c")
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "A.3" # ab x co
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "c" |
                                      genodat[idx[idx2],idx.odd] == "c" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "ac"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "c" |
                                      genodat[idx[idx2],idx.odd] == "c" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "bc"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "b"
      }
      
      idx2 <- which(genofile[idx,idxP2.1] == "a" & genofile[idx,idxP2.2] == "o" |
                      genofile[idx,idxP2.1] == "o" & genofile[idx,idxP2.2] == "a")
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "B1.5" # ab x ao
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "a" ,arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "b"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "b" |
                                      genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "ab"
      }
      
      idx2 <- which(genofile[idx,idxP2.1] == "a" & genofile[idx,idxP2.2] == "b" |
                      genofile[idx,idxP2.1] == "b" & genofile[idx,idxP2.2] == "a")
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "B3.7" # ab x ab
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "b" |
                                      genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "ab"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "a" ,arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "b" ,arr.ind = T)] <- "b"
      }
      
      idx2 <- which(genofile[idx,idxP2.1] == "c" & genofile[idx,idxP2.2] == "c")
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "D1.9" # ab x cc
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "c" |
                                      genodat[idx[idx2],idx.odd] == "c" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "ac"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "c" |
                                      genodat[idx[idx2],idx.odd] == "c" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "bc"
      }
      
      
      idx2 <- which(genofile[idx,idxP2.1] == "a" & genofile[idx,idxP2.2] == "a" )
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "D1.10" # ab x aa
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "b" |
                                      genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "ab"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "a" ,arr.ind = T)] <- "a"
      }
      
      idx2 <- which(genofile[idx,idxP2.1] == "o" & genofile[idx,idxP2.2] == "o")
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "D1.11" # ab x oo
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "b"
        
      }
    }
    idx <- which(genofile[,idxP1.1] == "a" & genofile[,idxP1.2] == "o" |
                   genofile[,idxP1.1] == "o" & genofile[,idxP1.2] == "a")
    
    if(length(idx)>0){
      idx2 <- which(genofile[idx,idxP2.1] == "b" & genofile[idx,idxP2.2] == "o" |
                      genofile[idx,idxP2.1] == "o" & genofile[idx,idxP2.2] == "b")
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "A.4" # ao x bo 
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "b"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "b" |
                                      genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "ab"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "o" ,arr.ind = T)] <- "o"
      }
      
      idx2 <- which(genofile[idx,idxP2.1] == "a" & genofile[idx,idxP2.2] == "b" |
                      genofile[idx,idxP2.1] == "b" & genofile[idx,idxP2.2] == "a")
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "B2.6" # ao x ab
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "a" ,arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "a" |
                                      genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "ab"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "b"
        
      }
      
      idx2 <- which(genofile[idx,idxP2.1] == "a" & genofile[idx,idxP2.2] == "o" |
                      genofile[idx,idxP2.1] == "o" & genofile[idx,idxP2.2] == "a")
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "C.8" # ao x ao
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "a" ,arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "a" |
                                      genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "o",arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "o",arr.ind = T)] <- "o"
      }
      
      idx2 <- which(genofile[idx,idxP2.1] == "o" & genofile[idx,idxP2.2] == "o" )
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "D1.13" # ao x oo
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "o" ,arr.ind = T)] <- "o"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "a"
      }
    }
    idx <- which(genofile[,idxP1.1] == "b" & genofile[,idxP1.2] == "o" |
                   genofile[,idxP1.1] == "o" & genofile[,idxP1.2] == "b")
    
    if(length(idx)>0){
      idx2 <- which(genofile[idx,idxP2.1] == "a" & genofile[idx,idxP2.2] == "a" )
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "D1.12" # bo x aa
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "a" |
                                      genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "ab"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "a"
      }
    }
    
    idx <- which(genofile[,idxP2.1] == "a" & genofile[,idxP2.2] == "b" |
                   genofile[,idxP2.1] == "b" & genofile[,idxP2.2] == "a")
    
    if(length(idx)>0){
      idx2 <- which(genofile[idx,idxP1.1] == "a" & genofile[idx,idxP1.2] == "a" )
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "D2.15" # aa x ab
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "a" ,arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "a" |
                                      genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "ab"
      }
      
      idx2 <- which(genofile[idx,idxP1.1] == "c" & genofile[idx,idxP1.2] == "c" )
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "D2.14" # cc x ab
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "c" &  genodat[idx[idx2],idx.even] == "a" |
                                      genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "c",arr.ind = T)] <- "ac"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "c" |
                                      genodat[idx[idx2],idx.odd] == "c" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "bc"
      }
      
      idx2 <- which(genofile[idx,idxP1.1] == "o" & genofile[idx,idxP1.2] == "o" )
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "D2.16" # oo x ab
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "b"
        
      }
    }
    
    idx <- which(genofile[,idxP2.1] == "o" & genofile[,idxP2.2] == "b" |
                   genofile[,idxP2.1] == "b" & genofile[,idxP2.2] == "o")
    
    if(length(idx)>0){
      idx2 <- which(genofile[idx,idxP1.1] == "a" & genofile[idx,idxP1.2] == "a" )
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "D2.17" # aa x bo
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "a" |
                                      genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "ab"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "a"
      }
    }
    
    
    idx <- which(genofile[,idxP2.1] == "a" & genofile[,idxP2.2] == "o" |
                   genofile[,idxP2.1] == "o" & genofile[,idxP2.2] == "a")
    
    if(length(idx)>0){
      idx2 <- which(genofile[idx,idxP1.1] == "o" & genofile[idx,idxP1.2] == "o" )
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "D2.18" # oo x ao
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "o" ,arr.ind = T)] <- "o"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "a"
      }
    }
  } else {
    idx <- which(genofile[,idxP1.1] == "a" & genofile[,idxP1.2] == "a")
    if(length(idx)>0){
      idx2 <- which(genofile[idx,idxP2.1] == "b" & genofile[idx,idxP2.2] == "b" )
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "A.H.B" # ab x ab
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "a" ,arr.ind = T)] <- "a"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "b" |
                                      genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "ab"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "b" ,arr.ind = T)] <- "b"
      }
      idx2 <- which(genofile[idx,idxP2.1] == "o" & genofile[idx,idxP2.2] == "o" )
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "C.A" # ao x ao
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "a" ,arr.ind = T)] <- "c"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "a" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "a",arr.ind = T)] <- "c"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "o" ,arr.ind = T)] <- "a"
      }
    }
    idx <- which(genofile[,idxP1.1] == "o" & genofile[,idxP1.2] == "o")
    if(length(idx)>0){
      idx2 <- which(genofile[idx,idxP2.1] == "b" & genofile[idx,idxP2.2] == "b" )
      if(length(idx2)>0) {
        segr.type.in[idx[idx2]] <- "D.B" # bo x bo
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "b" ,arr.ind = T)] <- "d"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "b" &  genodat[idx[idx2],idx.even] == "o" |
                                      genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "b",arr.ind = T)] <- "d"
        gt_matrix[idx[idx2],][which(genodat[idx[idx2],idx.odd] == "o" &  genodat[idx[idx2],idx.even] == "o" ,arr.ind = T)] <- "b"
      }
    }
  }
  
  pos.mis <- sample(1:length(gt_matrix),length(gt_matrix)*(miss.perc/100))
  gt_matrix[pos.mis] <- "-"
  
  head1 <- paste("data type", cross)
  n.ind <- dim(gt_matrix)[2]
  head2 <- paste(n.ind, n.mk, "0", "0", "0")
  ind.names <- unique(sapply(lapply(strsplit(colnames(genodat), split = "_"), "[",1:2), function(x) paste(x,collapse = "_")))
  mk.names <- genofile$marker
  geno <- cbind(paste0("*",mk.names), paste0(segr.type.in, "\t"), gt_matrix)
  
  write.table(c(head1, head2, paste(ind.names, collapse = " "), apply(geno, 1, function(x) paste(x, collapse = " "))), 
              file=out.file, row.names = F, col.names = F, quote = F)
  
}


