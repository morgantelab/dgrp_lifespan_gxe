###Load libraries needed
library(optparse)
library(BGLR)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--trm_file"), type="character")
parser <- add_option(parser, c("--pheno_file"), type="character")
parser <- add_option(parser, c("--cv_file"), type="character")
parser <- add_option(parser, c("--rep"), type="integer")
parser <- add_option(parser, c("--sex"), type="integer")
parser <- add_option(parser, c("--temp"), type="integer")
parser <- add_option(parser, c("--R2"), type="numeric")
parser <- add_option(parser, c("--niter"), type="integer")
parser <- add_option(parser, c("--nburn"), type="integer")
parser <- add_option(parser, c("--thin"), type="integer")
parser <- add_option(parser, c("--verbose"), type="logical")
parser <- add_option(parser, c("--output_file"), type="character")
parser <- add_option(parser, c("--intermediate_file"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

trm_file <- outparse$trm_file
pheno_file <- outparse$pheno_file
cv_file <- outparse$cv_file
repp <- outparse$rep
temp <- outparse$temp
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
TRM <- readRDS(trm_file)
dat <- readRDS(pheno_file)
cv <- readRDS(cv_file)

###Extract only temperature requested
obs_to_keep <- which(dat$temp == temp)
dat <- dat[obs_to_keep,]

###Order TRM according to phenotype data
TRM <- TRM[unique(dat$line_id), unique(dat$line_id)]

###Set up linear predictor for BGLR
ETA <- list(T=list(K=TRM,model="RKHS"))
  
###Set test set observations to missing
yNa <- dat$y
whichNa <- which(dat$line_id %in% cv[[repp]])
yNa[whichNa] <- NA

###Fit BGLR to the training set
fit <- BGLR(y=yNa, ETA=ETA, R2=R2, nIter=niter, burnIn=nburn, thin=thin, 
            saveAt=intermediate_file, verbose=FALSE)

###Predict phenotypes in the test set
yHat_test <- fit$yHat[whichNa]

###Compile results and write them to a file
res <- list(model_fit=fit, y_test=dat$y[whichNa], yhat_test=yHat_test)
saveRDS(res, output_file)



