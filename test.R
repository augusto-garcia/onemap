library(onemap)
library(reshape2)

df <- read_onemap("vcf_example_out.raw")

df <- create_probs(df, error = 10^-5, cross = "outcross")
str(df)

twopts <- rf_2pts(df)
seq1 <- make_seq(twopts,"all")
lgs <- group(seq1)
lg4 <- make_seq(lgs,4)

lg4.ord <- order_seq(lg4)
map.lg4.df <- make_seq(lg4.ord, "force")


df.f2 <- read_onemap("vcf_example_f2.raw")
df.f2 <- create_probs(df.f2, cross="f2")
str(df.f2)
head(df.f2$error)

twopts <- rf_2pts(df.f2)
seq1 <- make_seq(twopts,"all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)

lg2.ord <- order_seq(input.seq = lg2)
map.lg2.df <- make_seq(lg2.ord, "force")


df.bc <- read_onemap("vcf_example_bc.raw")
df.bc <- create_probs(df.bc, cross="backcross")
str(df.bc)
head(df.bc$error)

twopts <- rf_2pts(df.bc)
seq1 <- make_seq(twopts,"all")
lgs <- group(seq1)
lg2 <- make_seq(lgs,2)

lg2.ord <- order_seq(lg2)
map.lg2.df <- make_seq(lg2.ord, "force")

### Function
create_probs <- function(df, error = 10^(-5), cross=c("outcross", "backcross", "rils")){
  probs <- melt(t(df$geno))
  probs$type <- rep(df$segr.type.num, df$n.ind)
  
  # Global error according to observated data
  
  prob <- matrix(NA, nrow=length(probs$value), ncol = 4)
  
  idx <- which(probs$value == 0)
  prob[idx,] <- 1
  
  if(cross == "outcross"){
    # A
    idx <- which(probs$value == 1  & probs$type == 1)
    prob[idx,] <- c(rep(1- error, length(idx)), rep(error/3, 3*length(idx)))
    idx <- which(probs$value == 2  & probs$type == 1)
    prob[idx,] <- c(rep(error/3, length(idx)), rep(1-error, length(idx)), rep(error/3, 2*length(idx)))
    idx <- which(probs$value == 3  & probs$type == 1)
    prob[idx,] <- c(rep(error/3, 2*length(idx)), rep(1-error, length(idx)), rep(error/3, length(idx)))
    idx <- which(probs$value == 4  & probs$type == 1)
    prob[idx,] <- c(rep(error/3, 3*length(idx)), rep(1-error, length(idx)))
    
    # B1
    idx <- which(probs$value == 1  & probs$type == 2)
    prob[idx,] <- c(rep(1- error, 2*length(idx)), rep(error/2, 2*length(idx)))
    idx <- which(probs$value == 2  & probs$type == 2)
    prob[idx,] <- c(rep(error/3, 2*length(idx)), rep(1-error, length(idx)), rep(error/3, length(idx)))
    idx <- which(probs$value == 3  & probs$type == 2)
    prob[idx,] <- c(rep(error/3, 3*length(idx)), rep(1-error, length(idx)))
    
    # B2
    idx <- which(probs$value == 1  & probs$type == 3)
    prob[idx,] <- c(rep(1- error, length(idx)), rep(error/2, length(idx)), rep(1- error, length(idx)), rep(error/2, length(idx)))
    idx <- which(probs$value == 2  & probs$type == 3)
    prob[idx,] <- c(rep(error/3, length(idx)), rep(1-error, length(idx)), rep(error/3, 2*length(idx)))
    idx <- which(probs$value == 3  & probs$type == 3)
    prob[idx,] <- c(rep(error/3, 3*length(idx)), rep(1-error, length(idx)))
    
    # B3.7
    idx <- which(probs$value == 1  & probs$type == 4)
    prob[idx,] <- c(rep(1- error, length(idx)), rep(error/3, 3*length(idx)))
    idx <- which(probs$value == 2  & probs$type == 4)
    prob[idx,] <- c(rep(error/2, length(idx)), rep(1-error, 2*length(idx)), rep(error/2, length(idx)))
    idx <- which(probs$value == 3  & probs$type == 4)
    prob[idx,] <- c(rep(error/3, 3*length(idx)), rep(1-error, length(idx)))
    
    # C
    idx <- which(probs$value == 1  & probs$type == 5)
    prob[idx,] <- c(rep(1- error, 3*length(idx)), rep(error, length(idx)))
    idx <- which(probs$value == 2  & probs$type == 5)
    prob[idx,] <- c(rep(error/3, 3*length(idx)), rep(1-error, length(idx)))
    
    # D1
    idx <- which(probs$value == 1  & probs$type == 6)
    prob[idx,] <- c(rep(1- error, 2*length(idx)), rep(error/2, 2*length(idx)))
    idx <- which(probs$value == 2  & probs$type == 6)
    prob[idx,] <- c(rep(error/2, 2*length(idx)), rep(1-error, 2*length(idx)))
    
    # D2
    idx <- which(probs$value == 1  & probs$type == 7)
    prob[idx,] <- c(rep(1- error, length(idx)), rep(error/2, length(idx)), rep(1- error, length(idx)), rep(error/2, length(idx)))
    idx <- which(probs$value == 2  & probs$type == 7)
    prob[idx,] <- c(rep(error/2, length(idx)), rep(1-error, length(idx)), rep(error/2, length(idx)), rep(1-error, length(idx)))
    
  } else if(cross == "f2"){
    idx <- which(probs$value == 1)
    prob[idx,] <- c(rep(1- error, length(idx)), rep(error/3, 3*length(idx)))
    idx <- which(probs$value == 2)
    prob[idx,] <- c(rep(error/2, length(idx)), rep(1-error, 2*length(idx)), rep(error/2, length(idx)))
    idx <- which(probs$value == 3)
    prob[idx,] <- c(rep(error/3, 3*length(idx)), rep(1-error, length(idx)))
    idx <- which(probs$value == 4)
    prob[idx,] <- c(rep(1- error/3, 3*length(idx)), rep(error, length(idx)))
    idx <- which(probs$value == 5)
    prob[idx,] <- c(rep(error, length(idx)), rep(1-error/3, 3*length(idx)))
  } else if(cross == "backcross" | cross == "rils"){
    idx <- which(probs$value == 1)
    prob[idx,] <- c(rep(1- error, length(idx)), rep(error, length(idx)))
    idx <- which(probs$value == 2)
    prob[idx,] <- c(rep(error, length(idx)), rep(1-error, length(idx)))
    idx <- which(probs$value == 3)
    prob[idx,] <- c(rep(error, length(idx)), rep(1-error, length(idx)))
  }
  rownames(prob) <- paste0(probs$Var1, "_", probs$Var2)
  df$error <- prob
  return(df)
}
