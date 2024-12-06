###Load libraries needed
library(optparse)
library(dplyr)
library(BGLR)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--grm_file"), type="character")
parser <- add_option(parser, c("--pheno_file"), type="character")
parser <- add_option(parser, c("--cv_file"), type="character")
parser <- add_option(parser, c("--rep"), type="integer")
parser <- add_option(parser, c("--model"), type="character")
parser <- add_option(parser, c("--R2"), type="numeric")
parser <- add_option(parser, c("--niter"), type="integer")
parser <- add_option(parser, c("--nburn"), type="integer")
parser <- add_option(parser, c("--thin"), type="integer")
parser <- add_option(parser, c("--verbose"), type="logical")
parser <- add_option(parser, c("--output_file"), type="character")
parser <- add_option(parser, c("--intermediate_file"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

grm_file <- outparse$grm_file
pheno_file <- outparse$pheno_file
cv_file <- outparse$cv_file
repp <- outparse$rep
model <- outparse$model
R2 <- outparse$R2
niter <- outparse$niter
nburn <- outparse$nburn
thin <- outparse$thin
verbose <- outparse$verbose
output_file <- outparse$output_file
intermediate_file <- outparse$intermediate_file
seed <- outparse$seed

###Set seed
set.seed(seed)

###Load data
G <- readRDS(grm_file)
dat <- readRDS(pheno_file)
cv <- readRDS(cv_file)

###Create an indicator for the environments
env_groups <- dat %>% 
  group_by(temp, sex) %>%
  mutate(env = cur_group_id()) %>%
  ungroup %>%
  select(env) %>%
  as.data.frame
env_groups <- env_groups[,1]

###Order GRM according to phenotype data
G <- G[unique(dat$line_id), unique(dat$line_id)]

###Create ZGZ
if(model %in% c("G", "GandE", "GxE", "G_het_var", "GandE_het_var", "GxE_het_var")){
  Z <- model.matrix(y~line_id-1, dat)
  ZGZ <- Z%*%G%*%t(Z)
  
  if(model %in% c("G_het_var", "GandE_het_var", "GxE_het_var")){
    ZGZ_eig <- eigen(ZGZ)
    PC_ZGZ <- ZGZ_eig$vectors
    for(i in 1:ncol(PC_ZGZ)){  PC_ZGZ[,i] <- PC_ZGZ[,i]*sqrt(ZGZ_eig$values[i]) }
    PC_ZGZ <- PC_ZGZ[,ZGZ_eig$values>1e-5]
  }
}

###Create E
if(model %in% c("E", "GandE", "GxE", "E_het_var", "GandE_het_var", "GxE_het_var")){
  X <- dat[, c("temp", "sex")]
  We <- scale(X)
  E <- tcrossprod(We)
  E <- E/mean(diag(E))
  
  if(model %in% c("E_het_var", "GandE_het_var", "GxE_het_var")){
    E_eig <- eigen(E)
    PC_E <- E_eig$vectors
    for(i in 1:ncol(PC_E)){  PC_E[,i] <- PC_E[,i]*sqrt(E_eig$values[i]) }
    PC_E <- PC_E[,E_eig$values>1e-5]
  }
}

###Set up linear predictor for RKHS
if(model %in% c("G", "E", "GandE", "GxE")){
  if(model=="G"){
    
    ETA <- list(G=list(K=ZGZ,model="RKHS"))
    
  } else if(model=="E"){
    
    ETA <- list(E=list(K=E,model="RKHS"))
    
  } else if(model=="GandE"){
    
    ETA <- list(G=list(K=ZGZ,model="RKHS"),
                E=list(K=E,model="RKHS"))
    
  } else if(model=="GxE"){
    
    ZGZxE <- ZGZ*E
    ETA <- list(G=list(K=ZGZ,model="RKHS"),
                E=list(K=E,model="RKHS"),
                GxE=list(K=ZGZxE,model="RKHS"))
    
  }
  
  ###Set up linear predictor for BRR
} else if(model %in% c("G_het_var", "E_het_var", "GandE_het_var", "GxE_het_var")){
  if(model=="G_het_var"){
    
    ETA <- list(G=list(X=PC_ZGZ,model="BRR"))
    
  } else if(model=="E_het_var"){
    
    ETA <- list(E=list(X=PC_E,model="BRR"))
    
  } else if(model=="GandE_het_var"){
    
    ETA <- list(G=list(X=PC_ZGZ,model="BRR"),
                E=list(X=PC_E,model="BRR"))
    
  } else if(model=="GxE_het_var"){
    
    ZGZxE <- ZGZ*E
    
    ZGZxE_eig <- eigen(ZGZxE)
    PC_ZGZxE <- ZGZxE_eig$vectors
    for(i in 1:ncol(PC_ZGZxE)){  PC_ZGZxE[,i] <- PC_ZGZxE[,i]*sqrt(ZGZxE_eig$values[i]) }
    PC_ZGZxE <- PC_ZGZxE[,ZGZxE_eig$values>1e-5]
    
    ETA <- list(G=list(X=PC_ZGZ,model="BRR"),
                E=list(X=PC_E,model="BRR"),
                GxE=list(X=PC_ZGZxE,model="BRR"))
    
  }
}

###Set test set observations to missing
yNa <- dat$y
whichNa <- which(dat$obs_id %in% cv[[repp]])
yNa[whichNa] <- NA

###Fit BGLR to the training set
if(model %in% c("G", "E", "GandE", "GxE")){
  fit <- BGLR(y=yNa, ETA=ETA, R2=R2, nIter=niter, burnIn=nburn, thin=thin, 
              saveAt=intermediate_file, verbose=FALSE)
  
} else if(model %in% c("G_het_var", "E_het_var", "GandE_het_var", "GxE_het_var")){
  
  fit <- BGLR(y=yNa, ETA=ETA, R2=R2, nIter=niter, burnIn=nburn, thin=thin, groups=env_groups,
              saveAt=intermediate_file, verbose=FALSE)
}

###Predict phenotypes in the test set
yHat_test <- fit$yHat[whichNa]

###Compile results and write them to a file
res <- list(model_fit=fit, y_test=dat$y[whichNa], yhat_test=yHat_test)
saveRDS(res, output_file)



