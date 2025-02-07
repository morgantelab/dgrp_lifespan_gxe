###Load libraries needed
library(optparse)
library(genio)
library(dplyr)
library(BGLR)

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

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--grm_file"), type="character")
parser <- add_option(parser, c("--geno_file"), type="character")
parser <- add_option(parser, c("--pheno_file"), type="character")
parser <- add_option(parser, c("--cv_file"), type="character")
parser <- add_option(parser, c("--gwas_file"), type="character")
parser <- add_option(parser, c("--rep"), type="integer")
parser <- add_option(parser, c("--model"), type="character")
parser <- add_option(parser, c("--R2"), type="numeric")
parser <- add_option(parser, c("--niter"), type="integer")
parser <- add_option(parser, c("--nburn"), type="integer")
parser <- add_option(parser, c("--thin"), type="integer")
parser <- add_option(parser, c("--verbose"), type="logical")
parser <- add_option(parser, c("--gwas_p_thresh"), type="numeric")
parser <- add_option(parser, c("--output_file"), type="character")
parser <- add_option(parser, c("--intermediate_file"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

grm_file <- outparse$grm_file
geno_file <- outparse$geno_file
pheno_file <- outparse$pheno_file
cv_file <- outparse$cv_file
gwas_file <- outparse$gwas_file
repp <- outparse$rep
model <- outparse$model
R2 <- outparse$R2
niter <- outparse$niter
nburn <- outparse$nburn
thin <- outparse$thin
verbose <- outparse$verbose
gwas_p_thresh <- outparse$gwas_p_thresh
output_file <- outparse$output_file
intermediate_file <- outparse$intermediate_file
seed <- outparse$seed

###Set seed
set.seed(seed)

###Load data
G <- readRDS(grm_file)
geno <- read_plink(geno_file, verbose = FALSE)
dat <- readRDS(pheno_file)
cv <- readRDS(cv_file)
gwas <- readRDS(gwas_file)

###Create an indicator for the environments
env_groups <- dat %>% 
  group_by(temp, sex) %>%
  mutate(env = cur_group_id()) %>%
  ungroup %>%
  select(env) %>%
  as.data.frame
env_groups <- env_groups[,1]

###Order GRM and geno data according to phenotype data
G <- G[unique(dat$line_id), unique(dat$line_id)]
X <- geno$X[, unique(dat$line_id)]
rm(geno)

###Impute missing values and scale X
X <- mean_impute(X)
X <- colScale(X, center = TRUE, scale = TRUE, add_attr = FALSE)

###Select significant GxE interactions
to_sel <- which(gwas[, 9] < gwas_p_thresh)

if(length(to_sel) > 0){
  gwas_sel <- gwas[to_sel, ]
  
  ###Split X in significant and not significant variants for GxE
  X_sig_indexes <- which(colnames(X) %in% rownames(gwas_sel))
  X_sig <- X[, X_sig_indexes, drop=FALSE]
  X_not_sig <- X[, -X_sig_indexes, drop=FALSE]
  rm(X)
  
  ###Create two GxE GRMs
  G1 <- tcrossprod(X_sig)
  G1 <- G1/mean(diag(G1))
  
  G2 <- tcrossprod(X_not_sig)
  G2 <- G2/mean(diag(G2))
  
  ###Create ZGZ, ZG1Z and ZG2Z
  Z <- model.matrix(y~line_id-1, dat)
  ZGZ <- Z%*%G%*%t(Z)
  ZG1Z <- Z%*%G1%*%t(Z)
  ZG2Z <- Z%*%G2%*%t(Z)
  
  ###Create E
  Xe <- dat[, c("temp", "sex")]
  We <- scale(Xe)
  E <- tcrossprod(We)
  E <- E/mean(diag(E))
  
  ###Create ZG1ZxE and ZG2ZxE
  ZG1ZxE <- ZG1Z*E
  ZG2ZxE <- ZG2Z*E
  
  ###Set up linear predictor for RKHS
  if(model=="GxE_gwas"){
    
    ETA <- list(G=list(K=ZGZ,model="RKHS"),
                E=list(K=E,model="RKHS"),
                G1xE=list(K=ZG1ZxE,model="RKHS"),
                G2xE=list(K=ZG2ZxE,model="RKHS"))
    
    ###Set up linear predictor for BRR
  } else if(model=="GxE_gwas_het_var"){

    E_eig <- eigen(E)
    PC_E <- E_eig$vectors
    for(i in 1:ncol(PC_E)){  PC_E[,i] <- PC_E[,i]*sqrt(E_eig$values[i]) }
    PC_E <- PC_E[,E_eig$values>1e-5]
    
    ZGZ_eig <- eigen(ZGZ)
    PC_ZGZ <- ZGZ_eig$vectors
    for(i in 1:ncol(PC_ZGZ)){  PC_ZGZ[,i] <- PC_ZGZ[,i]*sqrt(ZGZ_eig$values[i]) }
    PC_ZGZ <- PC_ZGZ[,ZGZ_eig$values>1e-5]
    
    ZG1ZxE_eig <- eigen(ZG1ZxE)
    PC_ZG1ZxE <- ZG1ZxE_eig$vectors
    for(i in 1:ncol(PC_ZG1ZxE)){  PC_ZG1ZxE[,i] <- PC_ZG1ZxE[,i]*sqrt(ZG1ZxE_eig$values[i]) }
    PC_ZG1ZxE <- PC_ZG1ZxE[,ZG1ZxE_eig$values>1e-5]
    
    ZG2ZxE_eig <- eigen(ZG2ZxE)
    PC_ZG2ZxE <- ZG2ZxE_eig$vectors
    for(i in 1:ncol(PC_ZG2ZxE)){  PC_ZG2ZxE[,i] <- PC_ZG2ZxE[,i]*sqrt(ZG2ZxE_eig$values[i]) }
    PC_ZG2ZxE <- PC_ZG2ZxE[,ZG2ZxE_eig$values>1e-5]
    
    ETA <- list(G=list(X=PC_ZGZ,model="BRR"),
                E=list(X=PC_E,model="BRR"),
                G1xE=list(X=PC_ZG1ZxE,model="BRR"),
                G2xE=list(X=PC_ZG2ZxE,model="BRR"))
  }
  
  ###Set test set observations to missing
  yNa <- dat$y
  whichNa <- which(dat$obs_id %in% cv[[repp]])
  yNa[whichNa] <- NA
  
  gc()
  
  ###Fit BGLR to the training set
  if(model == "GxE_gwas"){
    fit <- BGLR(y=yNa, ETA=ETA, R2=R2, nIter=niter, burnIn=nburn, thin=thin, 
                saveAt=intermediate_file, verbose=FALSE)
    
  } else if(model == "GxE_gwas_het_var"){
    
    fit <- BGLR(y=yNa, ETA=ETA, R2=R2, nIter=niter, burnIn=nburn, thin=thin, groups=env_groups,
                saveAt=intermediate_file, verbose=FALSE)
  }
  
  ###Predict phenotypes in the test set
  yHat_test <- fit$yHat[whichNa]
  
  ###Compile results and write them to a file
  res <- list(model_fit=fit, y_test=dat$y[whichNa], yhat_test=yHat_test)
} else {
  res <- NULL
  warning("No GxE GWAS pvalues below the set threshold. Returning NULL object.")
}

saveRDS(res, output_file)


