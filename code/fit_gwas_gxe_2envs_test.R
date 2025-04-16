library(optparse)
library(genio)
library(rrBLUP)
library(dplyr)
library(pbapply)
library(parallel)

###Functions needed
mean_impute <- function(geno){
  f <- rowMeans(geno, na.rm = TRUE)
  for (i in 1:length(f)) {
    geno[i,][which(is.na(geno[i,]))] <- f[i]
  }
  return(geno)
}

colScale = function(x, center = TRUE, scale = TRUE, add_attr = TRUE, rows = NULL, cols = NULL) {
  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }
  ################
  # Get the column means
  ################
  cm = matrixStats::colMeans2(x, na.rm = TRUE)
  ################
  # Get the column sd
  ################
  if (scale) {
    csd = matrixStats::colSds(x, center = cm)
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x = t( (t(x) - cm) / csd )
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  return(x)
}

GxE_GWAS <- function(X, y, Z, G, E, C, mc.cores, method=c("grammar", "lmm", "lm")){
  
  ZX <- Z %*% X
  
  if(method=="grammar"){
    
    if(!is.null(C)){
      C <- cbind(int=1, c=C)
    } else {
      C <- matrix(1, nrow=length(y), ncol=1)
      colnames(C) <- "int"
    }
    
    fit_mix <- mixed.solve(y=y, X=cbind(C, E), Z=Z, K=G)
    res <- y - fit_mix$beta[1] - Z %*% fit_mix$u
    
    linreg <- function(i, X, y, E, C){
      
      gen <- rep(as.numeric(NA), 3)
      e1 <- rep(as.numeric(NA), 3)
      e2 <- rep(as.numeric(NA), 3)
      e1xe2 <- rep(as.numeric(NA), 3)
      gxe1 <- rep(as.numeric(NA), 3)
      gxe2 <- rep(as.numeric(NA), 3)
      gxe1xe2 <- rep(as.numeric(NA), 3)
      names(gen) <- names(e1) <- names(e2) <- names(e1xe2) <- names(gxe1) <- names(gxe2) <- names(gxe1xe2) <- c("b", "se", "pval")
      
      fit_lm <- lm(y ~ C + X[, i] + E[, 1] + E[, 2] + E[, 1]:E[, 2] + X[, i]:E[, 1] + X[, i]:E[, 2] + x:E[,1]:E[,2] - 1)
      fit_lm_sum <- summary(fit_lm)
      
      j <- ncol(C)
      
      gen[1] <- fit_lm_sum$coefficients[j+1, 1]
      gen[2] <- fit_lm_sum$coefficients[j+1, 2]
      gen[3] <- fit_lm_sum$coefficients[j+1, 4]
      e1[1] <- fit_lm_sum$coefficients[j+2, 1]
      e1[2] <- fit_lm_sum$coefficients[j+2, 2]
      e1[3] <- fit_lm_sum$coefficients[j+2, 4]
      e2[1] <- fit_lm_sum$coefficients[j+3, 1]
      e2[2] <- fit_lm_sum$coefficients[j+3, 2]
      e2[3] <- fit_lm_sum$coefficients[j+3, 4]
      e1xe1[1] <- fit_lm_sum$coefficients[j+4, 1]
      e1xe1[2] <- fit_lm_sum$coefficients[j+4, 2]
      e1xe1[3] <- fit_lm_sum$coefficients[j+4, 4]
      gxe1[1] <- fit_lm_sum$coefficients[j+5, 1]
      gxe1[2] <- fit_lm_sum$coefficients[j+5, 2]
      gxe1[3] <- fit_lm_sum$coefficients[j+5, 4]
      gxe2[1] <- fit_lm_sum$coefficients[j+6, 1]
      gxe2[2] <- fit_lm_sum$coefficients[j+6, 2]
      gxe2[3] <- fit_lm_sum$coefficients[j+6, 4]
      gxe1xe2[1] <- fit_lm_sum$coefficients[j+7, 1]
      gxe1xe2[2] <- fit_lm_sum$coefficients[j+7, 2]
      gxe1xe2[3] <- fit_lm_sum$coefficients[j+7, 4]
      
      return(list(g=gen, e1=e1, e2=e2, e1xe2=e1xe2, gxe1=gxe1, gxe2=gxe2, gxe1xe2=gxe1xe2))
    }
    
    p <- ncol(X)
    
    if(mc.cores>1){
      cl <- makeCluster(mc.cores)
      out <- parLapply(cl=cl, 1:p, linreg, ZX, res, E, C)
      stopCluster(cl)
    } else {
      out <- lapply(1:p, linreg, ZX, res, E, C)
    }
    
  } else if(method=="lmm"){
    
    lmmreg <- function(i, X, y, G, Z, E, C){
      
      gen <- rep(as.numeric(NA), 3)
      e1 <- rep(as.numeric(NA), 3)
      e2 <- rep(as.numeric(NA), 3)
      e1xe2 <- rep(as.numeric(NA), 3)
      gxe1 <- rep(as.numeric(NA), 3)
      gxe2 <- rep(as.numeric(NA), 3)
      gxe1xe2 <- rep(as.numeric(NA), 3)
      names(gen) <- names(e1) <- names(e2) <- names(e1xe2) <- names(gxe1) <- names(gxe2) <- names(gxe1xe2) <- c("b", "se", "pval")
      
      if(!is.null(C)){
        Xtest <- cbind(int=1, c=C, g=X[, i], e1=E[, 1], e2=E[, 2], e1xe2=E[, 1]*E[, 2], gxe1=X[, i] * E[, 1], 
                       gxe2=X[, i] * E[, 2],  gxe1xe2=X[, i] * E[, 1] * E[, 2])
      } else {
        Xtest <- cbind(int=1, g=X[, i], e1=E[, 1], e2=E[, 2], e1xe2=E[, 1]*E[, 2], gxe1=X[, i] * E[, 1], 
                       gxe2=X[, i] * E[, 2], gxe1xe2=X[, i] * E[, 1] * E[, 2])
      }
      
      fit_mix <- rrBLUP::mixed.solve(y=y, X=Xtest, Z=Z, K=G, SE=TRUE)
      
      j <- ncol(Xtest) - 7
      n <- nrow(Xtest)
      p <- ncol(Xtest)
      
      
      ###https://academic.oup.com/genetics/article/178/3/1709/6061473
      gen[1] <- fit_mix$beta[j+1]
      gen[2] <- fit_mix$beta.SE[j+1]
      gen[3] <- pf(q=(fit_mix$beta[j+1] / fit_mix$beta.SE[j+1])^2, df1 = 1, df2 = n-p, lower.tail = FALSE)
      e1[1] <- fit_mix$beta[j+2]
      e1[2] <- fit_mix$beta.SE[j+2]
      e1[3] <- pf(q=(fit_mix$beta[j+2] / fit_mix$beta.SE[j+2])^2, df1 = 1, df2 = n-p, lower.tail = FALSE)
      e2[1] <- fit_mix$beta[j+3]
      e2[2] <- fit_mix$beta.SE[j+3]
      e2[3] <- pf(q=(fit_mix$beta[j+3] / fit_mix$beta.SE[j+3])^2, df1 = 1, df2 = n-p, lower.tail = FALSE)
      e1xe2[1] <- fit_mix$beta[j+4]
      e1xe2[2] <- fit_mix$beta.SE[j+4]
      e1xe2[3] <- pf(q=(fit_mix$beta[j+4] / fit_mix$beta.SE[j+4])^2, df1 = 1, df2 = n-p, lower.tail = FALSE)
      gxe1[1] <- fit_mix$beta[j+5]
      gxe1[2] <- fit_mix$beta.SE[j+5]
      gxe1[3] <- pf(q=(fit_mix$beta[j+5] / fit_mix$beta.SE[j+5])^2, df1 = 1, df2 = n-p, lower.tail = FALSE)
      gxe2[1] <- fit_mix$beta[j+6]
      gxe2[2] <- fit_mix$beta.SE[j+6]
      gxe2[3] <- pf(q=(fit_mix$beta[j+6] / fit_mix$beta.SE[j+6])^2, df1 = 1, df2 = n-p, lower.tail = FALSE)
      gxe1xe2[1] <- fit_mix$beta[j+7]
      gxe1xe2[2] <- fit_mix$beta.SE[j+7]
      gxe1xe2[3] <- pf(q=(fit_mix$beta[j+7] / fit_mix$beta.SE[j+7])^2, df1 = 1, df2 = n-p, lower.tail = FALSE)
      
      
      return(list(g=gen, e1=e1, e2=e2, e1xe2=e1xe2, gxe1=gxe1, gxe2=gxe2, gxe1xe2=gxe1xe2))
    }
    
    p <- ncol(X)
    
    if(mc.cores>1){
      cl <- makeCluster(mc.cores)
      out <- parLapply(cl=cl, 1:p, lmmreg, ZX, y, G, Z, E, C)
      stopCluster(cl)
    } else {
      out <- lapply(1:p, lmmreg, ZX, y, G, Z, E, C)
    }
    
  } else if(method=="lm"){
    
    linreg <- function(i, X, y, E, C){
      
      gen <- rep(as.numeric(NA), 3)
      e1 <- rep(as.numeric(NA), 3)
      e2 <- rep(as.numeric(NA), 3)
      e1xe2 <- rep(as.numeric(NA), 3)
      gxe1 <- rep(as.numeric(NA), 3)
      gxe2 <- rep(as.numeric(NA), 3)
      gxe1xe2 <- rep(as.numeric(NA), 3)
      names(gen) <- names(e1) <- names(e2) <- names(e1xe2) <- names(gxe1) <- names(gxe2) <- names(gxe1xe2) <- c("b", "se", "pval")
      
      fit_lm <- lm(y ~ C + X[, i] + E[, 1] + E[, 2] + E[, 1]:E[, 2] + X[, i]:E[, 1] + X[, i]:E[, 2] + x:E[,1]:E[,2] - 1)
      fit_lm_sum <- summary(fit_lm)
      
      j <- ncol(C)
      
      gen[1] <- fit_lm_sum$coefficients[j+1, 1]
      gen[2] <- fit_lm_sum$coefficients[j+1, 2]
      gen[3] <- fit_lm_sum$coefficients[j+1, 4]
      e1[1] <- fit_lm_sum$coefficients[j+2, 1]
      e1[2] <- fit_lm_sum$coefficients[j+2, 2]
      e1[3] <- fit_lm_sum$coefficients[j+2, 4]
      e2[1] <- fit_lm_sum$coefficients[j+3, 1]
      e2[2] <- fit_lm_sum$coefficients[j+3, 2]
      e2[3] <- fit_lm_sum$coefficients[j+3, 4]
      e1xe1[1] <- fit_lm_sum$coefficients[j+4, 1]
      e1xe1[2] <- fit_lm_sum$coefficients[j+4, 2]
      e1xe1[3] <- fit_lm_sum$coefficients[j+4, 4]
      gxe1[1] <- fit_lm_sum$coefficients[j+5, 1]
      gxe1[2] <- fit_lm_sum$coefficients[j+5, 2]
      gxe1[3] <- fit_lm_sum$coefficients[j+5, 4]
      gxe2[1] <- fit_lm_sum$coefficients[j+6, 1]
      gxe2[2] <- fit_lm_sum$coefficients[j+6, 2]
      gxe2[3] <- fit_lm_sum$coefficients[j+6, 4]
      gxe1xe2[1] <- fit_lm_sum$coefficients[j+7, 1]
      gxe1xe2[2] <- fit_lm_sum$coefficients[j+7, 2]
      gxe1xe2[3] <- fit_lm_sum$coefficients[j+7, 4]
      
      return(list(g=gen, e1=e1, e2=e2, e1xe2=e1xe2, gxe1=gxe1, gxe2=gxe2, gxe1xe2=gxe1xe2))
    }
    
    p <- ncol(X)
    
    if(mc.cores>1){
      cl <- makeCluster(mc.cores)
      out <- parLapply(cl=cl, 1:p, linreg, ZX, y, E, C)
      stopCluster(cl)
    } else {
      out <- lapply(1:p, linreg, ZX, y, E, C)
    }
  }
  
  FUN <- function(x, i)
    do.call(rbind, lapply(x, `[[`, i))
  
  g_eff <- FUN(out, "g")
  colnames(g_eff) <- c("g_b", "g_se", "g_pval")
  e1_eff <- FUN(out, "e1")
  colnames(e1_eff) <- c("e1_b", "e1_se", "e1_pval")
  e2_eff <- FUN(out, "e2")
  colnames(e2_eff) <- c("e2_b", "e2_se", "e2_pval")
  e1xe2_eff <- FUN(out, "e1xe2")
  colnames(e1xe2_eff) <- c("e1xe2_b", "e1xe2_se", "e1xe2_pval")
  gxe1_eff <- FUN(out, "gxe1")
  colnames(gxe1_eff) <- c("gxe1_b", "gxe1_se", "gxe1_pval")
  gxe2_eff <- FUN(out, "gxe2")
  colnames(gxe2_eff) <- c("gxe2_b", "gxe2_se", "gxe2_pval")
  gxe1xe2_eff <- FUN(out, "gxe1xe2")
  colnames(gxe1xe2_eff) <- c("gxe1xe2_b", "gxe1xe2_se", "gxe1xe2_pval")
  
  dat_final <- cbind(g_eff, e1_eff, e2_eff, e1xe2_eff, gxe1_eff, gxe2_eff, gxe1xe2_eff)
  rownames(dat_final) <- colnames(X)
  
  return(dat_final)
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--geno_file"), type="character")
parser <- add_option(parser, c("--grm_file"), type="character")
parser <- add_option(parser, c("--pheno_file"), type="character")
parser <- add_option(parser, c("--output_file"), type="character")
parser <- add_option(parser, c("--method"), type="character")
parser <- add_option(parser, c("--standardize"), type="logical", default=FALSE)
parser <- add_option(parser, c("--ncores"), type="integer")
outparse <- parse_args(parser)

geno_file <- outparse$geno_file
grm_file <- outparse$grm_file
pheno_file <- outparse$pheno_file
output_file <- outparse$output_file
method <- outparse$method
ncores <- outparse$ncores
standardize <- outparse$standardize

###Set seed
set.seed(1)

###Load data
G <- readRDS(grm_file)
dat <- readRDS(pheno_file)
geno <- read_plink(geno_file, verbose = FALSE)

###Impute missing values
X <- mean_impute(geno$X)
X <- t(X)
rm(geno)

###Order GRM and geno data according to phenotype data
G <- G[unique(dat$line_id), unique(dat$line_id)]
X <- X[unique(dat$line_id), ]

###Create design matrix for random genetic effects
Z <- model.matrix(y~line_id-1, dat)

###Select environmental variables of interest
E <- cbind(sex=dat$sex, temp=dat$temp)

###Standardize variables, if requested
if(standardize){
  X <- colScale(X, center = TRUE, scale = TRUE, add_attr = FALSE)
  E <- colScale(E, center = TRUE, scale = TRUE, add_attr = FALSE)
}

###GxE GWAS
fit_gxe_gwas <- GxE_GWAS(X=X, y=dat$y, Z=Z, G=G, E=E, C=NULL, mc.cores=ncores, method=method)

###Compile results and write them to a file
saveRDS(fit_gxe_gwas, output_file)


