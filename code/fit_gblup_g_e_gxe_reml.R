###Load libraries needed
library(optparse)
library(dplyr)
library(sommer)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--grm_file"), type="character")
parser <- add_option(parser, c("--pheno_file"), type="character")
parser <- add_option(parser, c("--cv_file"), type="character")
parser <- add_option(parser, c("--rep"), type="integer")
parser <- add_option(parser, c("--model"), type="character")
parser <- add_option(parser, c("--verbose"), type="logical")
parser <- add_option(parser, c("--output_file"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

grm_file <- outparse$grm_file
pheno_file <- outparse$pheno_file
cv_file <- outparse$cv_file
repp <- outparse$rep
model <- outparse$model
verbose <- outparse$verbose
output_file <- outparse$output_file
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
  dplyr::select(env) %>%
  as.data.frame

dat$env_group <- env_groups[,1]


###Order GRM according to phenotype data
G <- G[unique(dat$line_id), unique(dat$line_id)]

###Create ZGZ
if(model %in% c("G", "GandE", "GxE", "G_het_var", "GandE_het_var", "GxE_het_var")){
  Z <- model.matrix(y~line_id-1, dat)
  rownames(Z) <- dat$obs_id
  ZGZ <- Z%*%G%*%t(Z)
}

###Create E
if(model %in% c("E", "GandE", "GxE", "E_het_var", "GandE_het_var", "GxE_het_var")){
  X <- dat[, c("temp", "sex")]
  We <- scale(X)
  E <- tcrossprod(We)
  E <- E/mean(diag(E))
  rownames(E) <- colnames(E) <- dat$obs_id
}

###Set kernels
if(model %in% c("G", "G_het_var")){
  
  rand <- "~ vsr(obs_id, Gu = ZGZ)"
  
} else if(model %in% c("E", "E_het_var")){
  
  rand <- "~vsr(obs_id, Gu = E)"
  
} else if(model %in% c("GandE", "GandE_het_var")){
  
  dat$obs_idd <- dat$obs_id
  
  rand <- "~ vsr(obs_id, Gu = ZGZ) + vsr(obs_idd, Gu = E)"
  
} else if(model %in% c("GxE", "GxE_het_var")){
  
  dat$obs_iddd <- dat$obs_idd <- dat$obs_id
  
  ZGZxE <- ZGZ*E
  
  rand <- "~ vsr(obs_id, Gu = ZGZ) + vsr(obs_idd, Gu = E) + vsr(obs_iddd, Gu = ZGZxE)"
}

###Set test set observations to missing
dat$yNa <- dat$y
whichNa <- which(dat$obs_id %in% cv[[repp]])
dat$yNa[whichNa] <- NA

###Fit REML to the training set
if(model %in% c("G", "E", "GandE", "GxE")){
  
  fit <- mmer(fixed=yNa ~ 1,
              random= as.formula(rand), 
              rcov= ~vsr(units),
              nIters=100,
              method="NR",
              dateWarning=FALSE, 
              date.warning=FALSE, 
              verbose=verbose,
              data=dat)
  
} else if(model %in% c("G_het_var", "E_het_var", "GandE_het_var", "GxE_het_var")){
  
  fit <- mmer(fixed=yNa ~ 1,
              random= as.formula(rand), 
              rcov= ~ vsr(dsr(as.character(env_group)), units),
              nIters=100,
              method="NR",
              dateWarning=FALSE, 
              date.warning=FALSE, 
              verbose=verbose,
              data=dat)
}  

###Predict phenotypes in the whole data and order them according to original order
U <- do.call("cbind", unlist(fit$U, recursive = FALSE))
yHat <- fit$Beta[1,3] + rowSums(U)
yHat <- yHat[dat$obs_id]

###Compile results and write them to a file
res <- list(model_fit=fit, y_test=dat$y[whichNa], yhat_test=yHat[whichNa])
saveRDS(res, output_file)
