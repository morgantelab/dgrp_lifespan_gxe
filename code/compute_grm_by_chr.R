###Load libraries needed
library(optparse)
library(genio)

###Functions needed
mean_impute <- function(geno){
  f <- rowMeans(geno, na.rm = TRUE)
  for (i in 1:length(f)) {
    geno[i,][which(is.na(geno[i,]))] <- f[i]
  }
  return(geno)
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--chr"), type="character")
outparse <- parse_args(parser)

geno_dat <- outparse$geno
output <- outparse$output
chr <- outparse$chr

set.seed(1)

###Read genotype data
geno <- read_plink(geno_dat, verbose = FALSE)

###Extract variants for selected chromosome
to_keep <- which(geno$bim$chr == chr)
  
###Impute missing values
X <- mean_impute(geno$X[to_keep, ])
X <- t(X)
  
###Compute GRM
W <- scale(X, center=TRUE, scale=TRUE)
G <- tcrossprod(W)
G <- G/mean(diag(G))
  
###Save GRM to file
saveRDS(G, file=output)


