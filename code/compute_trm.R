###Load libraries needed
library(optparse)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--expression"), type="character")
parser <- add_option(parser, c("--lines_to_keep"), type="character")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

exp_dat <- outparse$expression
lines_to_keep_dat <- outparse$lines_to_keep
output <- outparse$output

set.seed(1)

###Read data
exp <- read.table(exp_dat, header=TRUE, sep="\t", row.names = 1)
lines_to_keep <- read.table(lines_to_keep_dat, header=FALSE, sep="\t")

###Keep expression for only lines with phenotypes
exp <- exp[which(rownames(exp) %in% lines_to_keep[,1]),]

###Compute TRM
W <- scale(exp, center=TRUE, scale=TRUE)
TRM <- tcrossprod(W)
TRM <- TRM/mean(diag(TRM))

###Save GRM to file
saveRDS(TRM, file=output)
