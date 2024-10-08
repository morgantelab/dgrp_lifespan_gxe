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
}

saveRDS(test_sets, file=output)