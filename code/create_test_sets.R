###Load libraries needed
library(optparse)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--pheno"), type="character")
parser <- add_option(parser, c("--cv_scheme"), type="character")
parser <- add_option(parser, c("--n_reps"), type="integer")
parser <- add_option(parser, c("--test_set_prop"), type="numeric")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

pheno_file <- outparse$pheno
cv_scheme <- outparse$cv_scheme
n_reps <- outparse$n_reps
test_set_prop <- outparse$test_set_prop
output <- outparse$output
seed <- outparse$seed

###Set seed
set.seed(seed)

###Load data
dat <- readRDS(pheno_file)

###Assign lines to two groups
dgrp_lines <- 1:length(unique(dat$line_id))
lines_A <- sort(sample(x=dgrp_lines, size=length(dgrp_lines)*test_set_prop))
lines_B <- dgrp_lines[-lines_A]

lines_A <- dat$line_id[lines_A]
lines_B <- dat$line_id[lines_B]

###Create sets based on scenario
if(cv_scheme == "random_obs"){
  n <- nrow(dat)
  
  test_sets <- vector("list", n_reps)
  
  for(i in 1:n_reps){
    test_idx <- sort(sample(x=1:nrow(dat), size=round(n*test_set_prop)))
    test_sets[[i]] <- dat$obs_id[test_idx]
  }
} else if(cv_scheme == "random_lines"){
  lines_ids <- unique(dat$line_id)
  
  n <- length(lines_ids)
  test_sets <- vector("list", n_reps)
  
  for(i in 1:n_reps){
    test_lines_ids <- sample(x=lines_ids, size=round(n*test_set_prop))
    test_sets[[i]] <- dat$obs_id[which(dat$line_id %in% test_lines_ids)]
  }
  
} else if(cv_scheme == "random_within_env"){
  lines_ids <- unique(dat$line_id)
  
  n <- length(lines_ids)
  test_sets <- vector("list", n_reps)
  
  for(i in 1:n_reps){
    test_sets[[i]] <- sample(x=lines_ids, size=round(n*test_set_prop))
  }
} else if(cv_scheme == "temp"){
  test_sets <- vector("list", 3)
  
  it <- 0
  for(i in c(18,25,28)){
    it <- it+1
      
    test_idx <- which(dat$temp == i)
    test_sets[[it]] <- dat$obs_id[test_idx]
  }
} else if(cv_scheme == "temp_sex"){
  test_sets <- vector("list", 6)
  
  it <- 0
  for(i in c(18,25,28)){
    for(j in c(1,0)){
    it <- it+1
    
    test_idx <- which(dat$temp == i & dat$sex == j)
    test_sets[[it]] <- dat$obs_id[test_idx]
   }
  }
} else if(cv_scheme == "gradient_1_2"){
  test_sets <- vector("list", 2)
  
  test_sets[[1]] <- dat$obs_id[c(which((dat$line_id %in% lines_B) & dat$temp == 25), which(dat$temp == 28))]
  test_sets[[2]] <- dat$obs_id[c(which((dat$line_id %in% lines_A) & dat$temp == 25), which(dat$temp == 18))]
  
} else if(cv_scheme == "gradient_3_4"){
  test_sets <- vector("list", 2)
  
  test_sets[[1]] <- dat$obs_id[c(which((dat$sex == 1) & dat$temp == 25), which(dat$temp == 28))]
  test_sets[[2]] <- dat$obs_id[c(which((dat$se == 0) & dat$temp == 25), which(dat$temp == 18))]
  
}else if(cv_scheme == "sex"){
  test_sets <- vector("list", 2)
  
  test_sets[[1]] <- dat$obs_id[which(dat$sex == 1)]
  test_sets[[2]] <- dat$obs_id[which(dat$sex == 0)]
}

saveRDS(test_sets, file=output)