library(optparse)
library(genio)
library(rrBLUP)
library(dplyr)
library(parallel)

###Functions needed
mean_impute <- function(geno){
  f <- rowMeans(geno, na.rm = TRUE)
  for (i in 1:length(f)) {
    geno[i,][which(is.na(geno[i,]))] <- f[i]
  }
  return(t(geno))
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

GxE_GWAS <- function(X, y, Z, G, e, C, mc.cores, method=c("grammar", "lmm", "lm")){
  
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
      out <- parLapply(cl, 1:p, linreg, ZX, res, e, C)
      stopCluster(cl)
    } else {
      out <- lapply(1:p, linreg, ZX, res, e, C)
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
      
      gen[1] <- fit_mix$beta[j+1]
      gen[2] <- fit_mix$beta.SE[j+1]
      gen[3] <- pchisq(q=(fit_mix$beta[j+1] / fit_mix$beta.SE[j+1])^2, df=1, ncp = 0, lower.tail = FALSE)
      env[1] <- fit_mix$beta[j+2]
      env[2] <- fit_mix$beta.SE[j+2]
      env[3] <- pchisq(q=(fit_mix$beta[j+2] / fit_mix$beta.SE[j+2])^2, df=1, ncp = 0, lower.tail = FALSE)
      gxe[1] <- fit_mix$beta[j+3]
      gxe[2] <- fit_mix$beta.SE[j+3]
      gxe[3] <- pchisq(q=(fit_mix$beta[j+3] / fit_mix$beta.SE[j+3])^2, df=1, ncp = 0, lower.tail = FALSE)
      
      return(list(g=gen, e=env, gxe=gxe))
    }
    
    p <- ncol(X)
    
    if(mc.cores>1){
      cl <- makeCluster(mc.cores)
      out <- parLapply(cl, 1:p, lmmreg, ZX, y, G, Z, e, C)
      stopCluster(cl)
    } else {
      out <- lapply(1:p, lmmreg, ZX, y, G, Z, e, C)
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
      out <- parLapply(cl, 1:p, linreg, ZX, y, e, C)
      stopCluster(cl)
    } else {
      out <- lapply(1:p, linreg, ZX, y, e, C)
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
parser <- add_option(parser, c("--cv_file"), type="character")
parser <- add_option(parser, c("--output_file"), type="character")
parser <- add_option(parser, c("--method"), type="character")
parser <- add_option(parser, c("--environment"), type="character")
parser <- add_option(parser, c("--rep"), type="integer")
parser <- add_option(parser, c("--standardize"), type="logical", default=FALSE)
parser <- add_option(parser, c("--ncores"), type="integer")
outparse <- parse_args(parser)

geno_file <- outparse$geno_file
grm_file <- outparse$grm_file
pheno_file <- outparse$pheno_file
cv_file <- outparse$cv_file
output_file <- outparse$output_file
method <- outparse$method
envir <- outparse$environment
ncores <- outparse$ncores
standardize <- outparse$standardize
repp <- outparse$rep

###Set seed
set.seed(1)

###Load data
G <- readRDS(grm_file)
dat <- readRDS(pheno_file)
geno <- read_plink(geno_file, verbose = FALSE)
cv <- readRDS(cv_file)

###Compute environment means
env_means <- dat %>% group_by(sex, temp) %>% summarise(avg_y_by_env = mean(y)) %>% as.data.frame()
env_means <- env_means[order(env_means$avg_y_by_env),]
env_means$env_indicator <- 1:6

###Add environment means to main data
dat <- left_join(dat, env_means, by=join_by(sex,temp))

###Remove observations in the test
whichTest <- which(dat$obs_id %in% cv[[repp]])
dat <- dat[-whichTest, ]

###Order GRM and geno data according to phenotype data
G <- G[unique(dat$line_id), unique(dat$line_id)]
X <- geno$X[, unique(dat$line_id)]
rm(geno, env_means); gc()

###Impute missing values
X <- mean_impute(X)

###Create design matrix for random genetic effects
Z <- model.matrix(y~line_id-1, dat)

###Select environmental variable of interest
if(envir == "means"){
  e <- dat$avg_y_by_env
  C <- NULL
} else if(envir == "indicator"){
  e <- dat$env_indicator
  C <- NULL
} else if(envir == "sex"){
  e <- dat$sex
  C <- dat$temp
} else if(envir == "temp"){
  e <- dat$temp
  C <- dat$sex
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
fit_gxe_gwas <- GxE_GWAS(X=X, y=dat$y, Z=Z, G=G, e=e, C=C, mc.cores=ncores, method=method)

###Compile results and write them to a file
saveRDS(fit_gxe_gwas, output_file)
