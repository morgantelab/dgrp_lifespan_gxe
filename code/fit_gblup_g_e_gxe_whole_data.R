###Load libraries needed
library(optparse)
library(BGLR)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--grm_file"), type="character")
parser <- add_option(parser, c("--pheno_file"), type="character")
parser <- add_option(parser, c("--model"), type="character")
parser <- add_option(parser, c("--R2"), type="numeric")
parser <- add_option(parser, c("--niter"), type="integer")
parser <- add_option(parser, c("--nburn"), type="integer")
parser <- add_option(parser, c("--thin"), type="integer")
parser <- add_option(parser, c("--verbose"), type="logical", default=FALSE)
parser <- add_option(parser, c("--output_file"), type="character")
parser <- add_option(parser, c("--intermediate_file"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

grm_file <- outparse$grm_file
pheno_file <- outparse$pheno_file
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

###Order GRM according to phenotype data
G <- G[unique(dat$line_id), unique(dat$line_id)]

###Create ZGZ
if(model!="E"){
  Z <- model.matrix(y~line_id-1, dat)
  ZGZ <- Z%*%G%*%t(Z)
}

###Create E
if(model!="G"){
  X <- dat[, c("temp", "sex")]
  We <- scale(X)
  E <- tcrossprod(We)
  E <- E/mean(diag(E))
}

###Set up linear predictor for BGLR
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

###Fit BGLR to the training set
fit <- BGLR(y=dat$y, ETA=ETA, R2=R2, nIter=niter, burnIn=nburn, thin=thin, 
            saveAt=intermediate_file, verbose=verbose)

###Compile results and write them to a file
saveRDS(fit, output_file)



