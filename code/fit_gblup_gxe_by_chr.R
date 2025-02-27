###Load libraries needed
library(optparse)
library(dplyr)
library(BGLR)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--grm_file_1"), type="character")
parser <- add_option(parser, c("--grm_file_2"), type="character")
parser <- add_option(parser, c("--grm_file_3"), type="character")
parser <- add_option(parser, c("--grm_file_4"), type="character")
parser <- add_option(parser, c("--grm_file_5"), type="character")
parser <- add_option(parser, c("--grm_file_6"), type="character")
parser <- add_option(parser, c("--grm_file_all"), type="character")
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

grm_file_all <- outparse$grm_file_all
grm_file_1 <- outparse$grm_file_1
grm_file_2 <- outparse$grm_file_2
grm_file_3 <- outparse$grm_file_3
grm_file_4 <- outparse$grm_file_4
grm_file_5 <- outparse$grm_file_5
grm_file_6 <- outparse$grm_file_6
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
G_all <- readRDS(grm_file_all)
G_1 <- readRDS(grm_file_1)
G_2 <- readRDS(grm_file_2)
G_3 <- readRDS(grm_file_3)
G_4 <- readRDS(grm_file_4)
G_5 <- readRDS(grm_file_5)
G_6 <- readRDS(grm_file_6)
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
G_all <- G_all[unique(dat$line_id), unique(dat$line_id)]
G_1 <- G_1[unique(dat$line_id), unique(dat$line_id)]
G_2 <- G_2[unique(dat$line_id), unique(dat$line_id)]
G_3 <- G_3[unique(dat$line_id), unique(dat$line_id)]
G_4 <- G_4[unique(dat$line_id), unique(dat$line_id)]
G_5 <- G_5[unique(dat$line_id), unique(dat$line_id)]
G_6 <- G_6[unique(dat$line_id), unique(dat$line_id)]

###Create ZGZ
Z <- model.matrix(y~line_id-1, dat)
ZGZ_all <- Z%*%G_all%*%t(Z)
ZGZ_1 <- Z%*%G_1%*%t(Z)
ZGZ_2 <- Z%*%G_2%*%t(Z)
ZGZ_3 <- Z%*%G_3%*%t(Z)
ZGZ_4 <- Z%*%G_4%*%t(Z)
ZGZ_5 <- Z%*%G_5%*%t(Z)
ZGZ_6 <- Z%*%G_6%*%t(Z)
  
###Create E
X <- dat[, c("temp", "sex")]
We <- scale(X)
E <- tcrossprod(We)
E <- E/mean(diag(E))
  
###Create ZGZxE
ZGZ_1xE <- ZGZ_1*E
ZGZ_2xE <- ZGZ_2*E
ZGZ_3xE <- ZGZ_3*E
ZGZ_4xE <- ZGZ_4*E
ZGZ_5xE <- ZGZ_5*E
ZGZ_6xE <- ZGZ_6*E

if(model == "GxE_by_chr"){
  ETA <- list(G=list(K=ZGZ_all,model="RKHS"),
              E=list(K=E,model="RKHS"),
              G_1xE=list(K=ZGZ_1xE,model="RKHS"),
              G_2xE=list(K=ZGZ_2xE,model="RKHS"),
              G_3xE=list(K=ZGZ_3xE,model="RKHS"),
              G_4xE=list(K=ZGZ_4xE,model="RKHS"),
              G_5xE=list(K=ZGZ_5xE,model="RKHS"),
              G_6xE=list(K=ZGZ_6xE,model="RKHS"))
  
} else if(model == "GxE_by_chr_het_var"){
  
  ZGZ_all_eig <- eigen(ZGZ_all)
  PC_ZGZ_all <- ZGZ_all_eig$vectors
  for(i in 1:ncol(PC_ZGZ_all)){  PC_ZGZ_all[,i] <- PC_ZGZ_all[,i]*sqrt(ZGZ_all_eig$values[i]) }
  PC_ZGZ_all <- PC_ZGZ_all[,ZGZ_all_eig$values>1e-5]
  
  E_eig <- eigen(E)
  PC_E <- E_eig$vectors
  for(i in 1:ncol(PC_E)){  PC_E[,i] <- PC_E[,i]*sqrt(E_eig$values[i]) }
  PC_E <- PC_E[,E_eig$values>1e-5]
  
  ZGZ_1xE_eig <- eigen(ZGZ_1xE)
  PC_ZGZ_1xE <- ZGZ_1xE_eig$vectors
  for(i in 1:ncol(PC_ZGZ_1xE)){  PC_ZGZ_1xE[,i] <- PC_ZGZ_1xE[,i]*sqrt(ZGZ_1xE_eig$values[i]) }
  PC_ZGZ_1xE <- PC_ZGZ_1xE[,ZGZ_1xE_eig$values>1e-5]
  
  ZGZ_2xE_eig <- eigen(ZGZ_2xE)
  PC_ZGZ_2xE <- ZGZ_2xE_eig$vectors
  for(i in 1:ncol(PC_ZGZ_2xE)){  PC_ZGZ_2xE[,i] <- PC_ZGZ_2xE[,i]*sqrt(ZGZ_2xE_eig$values[i]) }
  PC_ZGZ_2xE <- PC_ZGZ_2xE[,ZGZ_2xE_eig$values>1e-5]
  
  ZGZ_3xE_eig <- eigen(ZGZ_3xE)
  PC_ZGZ_3xE <- ZGZ_3xE_eig$vectors
  for(i in 1:ncol(PC_ZGZ_3xE)){  PC_ZGZ_3xE[,i] <- PC_ZGZ_3xE[,i]*sqrt(ZGZ_3xE_eig$values[i]) }
  PC_ZGZ_3xE <- PC_ZGZ_3xE[,ZGZ_3xE_eig$values>1e-5]
  
  ZGZ_4xE_eig <- eigen(ZGZ_4xE)
  PC_ZGZ_4xE <- ZGZ_4xE_eig$vectors
  for(i in 1:ncol(PC_ZGZ_4xE)){  PC_ZGZ_4xE[,i] <- PC_ZGZ_4xE[,i]*sqrt(ZGZ_4xE_eig$values[i]) }
  PC_ZGZ_4xE <- PC_ZGZ_4xE[,ZGZ_4xE_eig$values>1e-5]
  
  ZGZ_5xE_eig <- eigen(ZGZ_5xE)
  PC_ZGZ_5xE <- ZGZ_5xE_eig$vectors
  for(i in 1:ncol(PC_ZGZ_5xE)){  PC_ZGZ_5xE[,i] <- PC_ZGZ_5xE[,i]*sqrt(ZGZ_5xE_eig$values[i]) }
  PC_ZGZ_5xE <- PC_ZGZ_5xE[,ZGZ_5xE_eig$values>1e-5]
  
  ZGZ_6xE_eig <- eigen(ZGZ_6xE)
  PC_ZGZ_6xE <- ZGZ_6xE_eig$vectors
  for(i in 1:ncol(PC_ZGZ_6xE)){  PC_ZGZ_6xE[,i] <- PC_ZGZ_6xE[,i]*sqrt(ZGZ_6xE_eig$values[i]) }
  PC_ZGZ_6xE <- PC_ZGZ_6xE[,ZGZ_6xE_eig$values>1e-5]
  
  ETA <- list(G=list(X=PC_ZGZ_all,model="BRR"),
              E=list(X=PC_E,model="BRR"),
              G_1xE=list(X=PC_ZGZ_1xE,model="BRR"),
              G_2xE=list(X=PC_ZGZ_2xE,model="BRR"),
              G_3xE=list(X=PC_ZGZ_3xE,model="BRR"),
              G_4xE=list(X=PC_ZGZ_4xE,model="BRR"),
              G_5xE=list(X=PC_ZGZ_5xE,model="BRR"),
              G_6xE=list(X=PC_ZGZ_6xE,model="BRR"))
}


###Set test set observations to missing
yNa <- dat$y
whichNa <- which(dat$obs_id %in% cv[[repp]])
yNa[whichNa] <- NA

###Fit BGLR to the training set
if(model == "GxE_by_chr"){
  fit <- BGLR(y=yNa, ETA=ETA, R2=R2, nIter=niter, burnIn=nburn, thin=thin, 
              saveAt=intermediate_file, verbose=FALSE)
  
} else if(model == "GxE_by_chr_het_var"){
  
  fit <- BGLR(y=yNa, ETA=ETA, R2=R2, nIter=niter, burnIn=nburn, thin=thin, groups=env_groups,
              saveAt=intermediate_file, verbose=FALSE)
}

###Predict phenotypes in the test set
yHat_test <- fit$yHat[whichNa]

###Compile results and write them to a file
res <- list(model_fit=fit, y_test=dat$y[whichNa], yhat_test=yHat_test)
saveRDS(res, output_file)



