###Load libraries needed
library(optparse)
library(tidyr)
library(BGLR)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--grm_file"), type="character")
parser <- add_option(parser, c("--pheno_file"), type="character")
parser <- add_option(parser, c("--R2"), type="numeric")
parser <- add_option(parser, c("--ResCov_type"), type="character")
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
R2 <- outparse$R2
niter <- outparse$niter
nburn <- outparse$nburn
thin <- outparse$thin
verbose <- outparse$verbose
output_file <- outparse$output_file
intermediate_file <- outparse$intermediate_file
ResCov_type <- outparse$ResCov_type ##c("UN","DIAG","REC","FA")
seed <- outparse$seed

###Set seed
set.seed(seed)

###Load data
G <- readRDS(grm_file)
dat <- readRDS(pheno_file)

dat_wide <- dat %>%
  pivot_wider(names_from = c(temp, sex), values_from = c(y,obs_id)) %>% as.data.frame()
rownames(dat_wide) <- dat_wide$line_id

###Order GRM according to phenotype data
G <- G[dat_wide$line_id, dat_wide$line_id]

###Extract phenotype matrix
y <- as.matrix(dat_wide[, 2:7])

###Fit BGLR to the training set
ETA <- list(G=list(K=G,model="RKHS"))

fit <- Multitrait(y=y, ETA=ETA, R2=R2, nIter=niter, burnIn=nburn, thin=thin, 
                  resCov=list(type=ResCov_type), saveAt=intermediate_file, verbose=verbose)

###Compile results and write them to a file
saveRDS(fit, output_file)



