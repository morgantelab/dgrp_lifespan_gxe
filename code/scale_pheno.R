###Load library
library(optparse)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--scale_y"), type="character")
parser <- add_option(parser, c("--input"), type="character")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

scale_y <- outparse$scale_y
input <- outparse$input
output <- outparse$output

###Set seed
set.seed(1)

###Read data
dat <- readRDS(input)

if(scale_y == "scaled_within_env"){
  library(dplyr)
  
  ###Scale phenotype within each sex/temp combo --> variance is very different
  dat_scaled <- dat %>% group_by(sex, temp) %>% mutate(y_s = scale(y, center=TRUE, scale=TRUE)) %>% as.data.frame()
  dat$y <- dat_scaled$y_s
} else if(scale_y == "scaled_overall"){
  dat$y <- scale(dat$y, center=TRUE, scale=TRUE)
}

###Save data
saveRDS(dat, file=output)

