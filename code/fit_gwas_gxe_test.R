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

GxE_GWAS_F <- function(X, y, Z, G, e, C, mc.cores, method=c("grammar", "lmm", "lm")){
  
  ZX <- Z %*% X
  
  if(method=="grammar"){
    
    if(!is.null(C)){
      C <- cbind(int=1, c=C)
    } else {
      C <- matrix(1, nrow=length(y), ncol=1)
      colnames(C) <- "int"
    }
    
    fit_mix <- mixed.solve(y=y, X=cbind(C, e), Z=Z, K=G)
    res <- y - fit_mix$beta[1] - Z %*% fit_mix$u
    
    linreg <- function(i, X, y, e, C){
      
      gen <- rep(as.numeric(NA), 3)
      env <- rep(as.numeric(NA), 3)
      gxe <- rep(as.numeric(NA), 3)
      names(gen) <- names(env) <- names(gxe) <- c("b", "se", "pval")
      
      fit_lm <- lm(y ~ C + X[, i]*e - 1)
      fit_lm_sum <- summary(fit_lm)
      
      j <- ncol(C)
      
      gen[1] <- fit_lm_sum$coefficients[j+1, 1]
      gen[2] <- fit_lm_sum$coefficients[j+1, 2]
      gen[3] <- fit_lm_sum$coefficients[j+1, 4]
      env[1] <- fit_lm_sum$coefficients[j+2, 1]
      env[2] <- fit_lm_sum$coefficients[j+2, 2]
      env[3] <- fit_lm_sum$coefficients[j+2, 4]
      gxe[1] <- fit_lm_sum$coefficients[j+3, 1]
      gxe[2] <- fit_lm_sum$coefficients[j+3, 2]
      gxe[3] <- fit_lm_sum$coefficients[j+3, 4]
      
      return(list(g=gen, e=env, gxe=gxe))
    }
    
    p <- ncol(X)
    
    if(mc.cores>1){
      cl <- makeCluster(mc.cores)
      out <- pblapply(cl=cl, 1:p, linreg, ZX, res, e, C)
      stopCluster(cl)
    } else {
      out <- pblapply(1:p, linreg, ZX, res, e, C)
    }
  
  } else if(method=="lmm"){
    
    lmmreg <- function(i, X, y, G, Z, e, C){
      
      gen <- rep(as.numeric(NA), 3)
      env <- rep(as.numeric(NA), 3)
      gxe <- rep(as.numeric(NA), 3)
      names(gen) <- names(env) <- names(gxe) <- c("b", "se", "pval")
      
      if(!is.null(C)){
        Xtest <- cbind(int=1, c=C, g=X[, i], e=e, gxe=X[, i] * e)
      } else {
        Xtest <- cbind(int=1, g=X[, i], e=e, gxe=X[, i] * e)
      }
      
      fit_mix <- rrBLUP::mixed.solve(y=y, X=Xtest, Z=Z, K=G, SE=TRUE)
      
      j <- ncol(Xtest) - 3
      n <- nrow(Xtest)
      p <- ncol(Xtest)
      
      
      ###https://academic.oup.com/genetics/article/178/3/1709/6061473
      gen[1] <- fit_mix$beta[j+1]
      gen[2] <- fit_mix$beta.SE[j+1]
      gen[3] <- pf(q=(fit_mix$beta[j+1] / fit_mix$beta.SE[j+1])^2, df1 = 1, df2 = n-p, lower.tail = FALSE)
      env[1] <- fit_mix$beta[j+2]
      env[2] <- fit_mix$beta.SE[j+2]
      env[3] <- pf(q=(fit_mix$beta[j+2] / fit_mix$beta.SE[j+2])^2, df1 = 1, df2 = n-p, lower.tail = FALSE)
      gxe[1] <- fit_mix$beta[j+3]
      gxe[2] <- fit_mix$beta.SE[j+3]
      gxe[3] <- pf(q=(fit_mix$beta[j+3] / fit_mix$beta.SE[j+3])^2, df1 = 1, df2 = n-p, lower.tail = FALSE)
      
      return(list(g=gen, e=env, gxe=gxe))
    }
    
    p <- ncol(X)
    
    if(mc.cores>1){
      cl <- makeCluster(mc.cores)
      out <- pblapply(cl=cl, 1:p, lmmreg, ZX, y, G, Z, e, C)
      stopCluster(cl)
    } else {
      out <- pblapply(1:p, lmmreg, ZX, y, G, Z, e, C)
    }
    
  } else if(method=="lm"){
    
    linreg <- function(i, X, y, e, C){
      
      gen <- rep(as.numeric(NA), 3)
      env <- rep(as.numeric(NA), 3)
      gxe <- rep(as.numeric(NA), 3)
      names(gen) <- names(env) <- names(gxe) <- c("b", "se", "pval")
      
      if(!is.null(C)){
        C <- cbind(int=1, c=C)
      } else {
        C <- matrix(1, nrow=length(y), ncol=1)
        colnames(C) <- "int"
      }
      
      fit_lm <- lm(y ~ C + X[, i]*e - 1)
      fit_lm_sum <- summary(fit_lm)
      
      j <- ncol(C)
      
      gen[1] <- fit_lm_sum$coefficients[j+1, 1]
      gen[2] <- fit_lm_sum$coefficients[j+1, 2]
      gen[3] <- fit_lm_sum$coefficients[j+1, 4]
      env[1] <- fit_lm_sum$coefficients[j+2, 1]
      env[2] <- fit_lm_sum$coefficients[j+2, 2]
      env[3] <- fit_lm_sum$coefficients[j+2, 4]
      gxe[1] <- fit_lm_sum$coefficients[j+3, 1]
      gxe[2] <- fit_lm_sum$coefficients[j+3, 2]
      gxe[3] <- fit_lm_sum$coefficients[j+3, 4]
      
      return(list(g=gen, e=env, gxe=gxe))
    }
    
    p <- ncol(X)
    
    if(mc.cores>1){
      cl <- makeCluster(mc.cores)
      out <- pblapply(cl=cl, 1:p, linreg, ZX, y, e, C)
      stopCluster(cl)
    } else {
      out <- pblapply(1:p, linreg, ZX, y, e, C)
    }
  }
  
  FUN <- function(x, i)
    do.call(rbind, lapply(x, `[[`, i))
  
  g_eff <- FUN(out, "g")
  colnames(g_eff) <- c("g_b", "g_se", "g_pval")
  e_eff <- FUN(out, "e")
  colnames(e_eff) <- c("e_b", "e_se", "e_pval")
  gxe_eff <- FUN(out, "gxe")
  colnames(gxe_eff) <- c("gxe_b", "gxe_se", "gxe_pval")
  
  dat_final <- cbind(g_eff, e_eff, gxe_eff)
  rownames(dat_final) <- colnames(X)
  
  return(dat_final)
}

GxE_GWAS_Chi2 <- function(X, y, Z, G, e, C, mc.cores, method=c("grammar", "lmm", "lm")){
  
  ZX <- Z %*% X
  
  if(method=="grammar"){
    
    if(!is.null(C)){
      C <- cbind(int=1, c=C)
    } else {
      C <- matrix(1, nrow=length(y), ncol=1)
      colnames(C) <- "int"
    }
    
    fit_mix <- mixed.solve(y=y, X=cbind(C, e), Z=Z, K=G)
    res <- y - fit_mix$beta[1] - Z %*% fit_mix$u
    
    linreg <- function(i, X, y, e, C){
      
      gen <- rep(as.numeric(NA), 3)
      env <- rep(as.numeric(NA), 3)
      gxe <- rep(as.numeric(NA), 3)
      names(gen) <- names(env) <- names(gxe) <- c("b", "se", "pval")
      
      fit_lm <- lm(y ~ C + X[, i]*e - 1)
      fit_lm_sum <- summary(fit_lm)
      
      j <- ncol(C)
      
      gen[1] <- fit_lm_sum$coefficients[j+1, 1]
      gen[2] <- fit_lm_sum$coefficients[j+1, 2]
      gen[3] <- fit_lm_sum$coefficients[j+1, 4]
      env[1] <- fit_lm_sum$coefficients[j+2, 1]
      env[2] <- fit_lm_sum$coefficients[j+2, 2]
      env[3] <- fit_lm_sum$coefficients[j+2, 4]
      gxe[1] <- fit_lm_sum$coefficients[j+3, 1]
      gxe[2] <- fit_lm_sum$coefficients[j+3, 2]
      gxe[3] <- fit_lm_sum$coefficients[j+3, 4]
      
      return(list(g=gen, e=env, gxe=gxe))
    }
    
    p <- ncol(X)
    
    if(mc.cores>1){
      cl <- makeCluster(mc.cores)
      out <- pblapply(cl=cl, 1:p, linreg, ZX, res, e, C)
      stopCluster(cl)
    } else {
      out <- pblapply(1:p, linreg, ZX, res, e, C)
    }
    
  } else if(method=="lmm"){
    
    lmmreg <- function(i, X, y, G, Z, e, C){
      
      gen <- rep(as.numeric(NA), 3)
      env <- rep(as.numeric(NA), 3)
      gxe <- rep(as.numeric(NA), 3)
      names(gen) <- names(env) <- names(gxe) <- c("b", "se", "pval")
      
      if(!is.null(C)){
        Xtest <- cbind(int=1, c=C, g=X[, i], e=e, gxe=X[, i] * e)
      } else {
        Xtest <- cbind(int=1, g=X[, i], e=e, gxe=X[, i] * e)
      }
      
      fit_mix <- rrBLUP::mixed.solve(y=y, X=Xtest, Z=Z, K=G, SE=TRUE)
      
      j <- ncol(Xtest) - 3
      n <- nrow(Xtest)
      p <- ncol(Xtest)
      
      gen[1] <- fit_mix$beta[j+1]
      gen[2] <- fit_mix$beta.SE[j+1]
      gen[3] <- pchisq(q=(fit_mix$beta[j+1] / fit_mix$beta.SE[j+1])^2, df = 1, lower.tail = FALSE)
      env[1] <- fit_mix$beta[j+2]
      env[2] <- fit_mix$beta.SE[j+2]
      env[3] <- pchisq(q=(fit_mix$beta[j+2] / fit_mix$beta.SE[j+2])^2, df = 1, lower.tail = FALSE)
      gxe[1] <- fit_mix$beta[j+3]
      gxe[2] <- fit_mix$beta.SE[j+3]
      gxe[3] <- pchisq(q=(fit_mix$beta[j+3] / fit_mix$beta.SE[j+3])^2, df = 1, lower.tail = FALSE)
      
      return(list(g=gen, e=env, gxe=gxe))
    }
    
    p <- ncol(X)
    
    if(mc.cores>1){
      cl <- makeCluster(mc.cores)
      out <- pblapply(cl=cl, 1:p, lmmreg, ZX, y, G, Z, e, C)
      stopCluster(cl)
    } else {
      out <- pblapply(1:p, lmmreg, ZX, y, G, Z, e, C)
    }
    
  } else if(method=="lm"){
    
    linreg <- function(i, X, y, e, C){
      
      gen <- rep(as.numeric(NA), 3)
      env <- rep(as.numeric(NA), 3)
      gxe <- rep(as.numeric(NA), 3)
      names(gen) <- names(env) <- names(gxe) <- c("b", "se", "pval")
      
      if(!is.null(C)){
        C <- cbind(int=1, c=C)
      } else {
        C <- matrix(1, nrow=length(y), ncol=1)
        colnames(C) <- "int"
      }
      
      fit_lm <- lm(y ~ C + X[, i]*e - 1)
      fit_lm_sum <- summary(fit_lm)
      
      j <- ncol(C)
      
      gen[1] <- fit_lm_sum$coefficients[j+1, 1]
      gen[2] <- fit_lm_sum$coefficients[j+1, 2]
      gen[3] <- fit_lm_sum$coefficients[j+1, 4]
      env[1] <- fit_lm_sum$coefficients[j+2, 1]
      env[2] <- fit_lm_sum$coefficients[j+2, 2]
      env[3] <- fit_lm_sum$coefficients[j+2, 4]
      gxe[1] <- fit_lm_sum$coefficients[j+3, 1]
      gxe[2] <- fit_lm_sum$coefficients[j+3, 2]
      gxe[3] <- fit_lm_sum$coefficients[j+3, 4]
      
      return(list(g=gen, e=env, gxe=gxe))
    }
    
    p <- ncol(X)
    
    if(mc.cores>1){
      cl <- makeCluster(mc.cores)
      out <- pblapply(cl=cl, 1:p, linreg, ZX, y, e, C)
      stopCluster(cl)
    } else {
      out <- pblapply(1:p, linreg, ZX, y, e, C)
    }
  }
  
  FUN <- function(x, i)
    do.call(rbind, lapply(x, `[[`, i))
  
  g_eff <- FUN(out, "g")
  colnames(g_eff) <- c("g_b", "g_se", "g_pval")
  e_eff <- FUN(out, "e")
  colnames(e_eff) <- c("e_b", "e_se", "e_pval")
  gxe_eff <- FUN(out, "gxe")
  colnames(gxe_eff) <- c("gxe_b", "gxe_se", "gxe_pval")
  
  dat_final <- cbind(g_eff, e_eff, gxe_eff)
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
parser <- add_option(parser, c("--environment"), type="character")
parser <- add_option(parser, c("--standardize"), type="logical", default=FALSE)
parser <- add_option(parser, c("--ncores"), type="integer")
outparse <- parse_args(parser)

geno_file <- outparse$geno_file
grm_file <- outparse$grm_file
pheno_file <- outparse$pheno_file
output_file <- outparse$output_file
method <- outparse$method
envir <- outparse$environment
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

###Compute environment means
env_means <- dat %>% group_by(sex, temp) %>% summarise(avg_y_by_env = mean(y)) %>% as.data.frame()
env_means <- env_means[order(env_means$avg_y_by_env),]
env_means$env_indicator <- 1:6

###Add environment means to main data
dat_final <- left_join(dat, env_means, by=join_by(sex,temp))
rm(dat, env_means); gc()

###Select environmental variable of interest
if(envir == "means"){
  e <- dat_final$avg_y_by_env
  C <- NULL
} else if(envir == "indicator"){
  e <- dat_final$env_indicator
  C <- NULL
} else if(envir == "sex"){
  e <- dat_final$sex
  C <- dat_final$temp
} else if(envir == "temp"){
  e <- dat_final$temp
  C <- dat_final$sex
}

###Standardize variables, if requested
if(standardize){
  X <- colScale(X, center = TRUE, scale = TRUE, add_attr = FALSE)
  e <- (e - mean(e, na.rm=TRUE))/sd(e, na.rm=TRUE)
  if(!is.null(C)){
    C <- (C - mean(C, na.rm=TRUE))/sd(C, na.rm=TRUE)
  }
}

###GxE GWAS
fit_gxe_gwas_F <- GxE_GWAS_F(X=X[,1:10], y=dat_final$y, Z=Z, G=G, e=e, C=C, mc.cores=ncores, method=method)

fit_gxe_gwas_Chi2 <- GxE_GWAS_Chi2(X=X[,1:10], y=dat_final$y, Z=Z, G=G, e=e, C=C, mc.cores=ncores, method=method)

###Compile results and write them to a file
saveRDS(fit_gxe_gwas, output_file)
 