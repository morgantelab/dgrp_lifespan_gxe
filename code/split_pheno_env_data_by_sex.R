###Load libraries needed
library(optparse)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--pheno"), type="character")
parser <- add_option(parser, c("--trm"), type="character")
parser <- add_option(parser, c("--sex"), type="integer")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

pheno_file <- outparse$pheno
trm_file <- outparse$trm
sex <- outparse$sex
output <- outparse$output

###Set seed
set.seed(1)

###Load data
TRM <- readRDS(trm_file)
dat <- readRDS(pheno_file)

###Filter data 
to_keep <- rownames(TRM)

dat <- dat[which(dat$line_id %in% to_keep), ]
dat <- dat[which(dat$sex %in% sex), ]

##Save output
saveRDS(dat, file=output)