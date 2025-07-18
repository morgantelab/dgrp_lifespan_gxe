###Load library
library(optparse)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--output_data"), type="character")
parser <- add_option(parser, c("--output_lines_ids"), type="character")
outparse <- parse_args(parser)

output_data <- outparse$output_data
output_lines_ids <- outparse$output_lines_ids

###Set seed
set.seed(1)

###Read in data
dat_F_18 <- read.table("/data2/morgante_lab/data/dgrp/phenotypes/lifespan/female_18c_mean_pheno.txt", header=FALSE, sep="\t")
dat_F_25 <- read.table("/data2/morgante_lab/data/dgrp/phenotypes/lifespan/female_25c_mean_pheno.txt", header=FALSE, sep="\t")
dat_F_28 <- read.table("/data2/morgante_lab/data/dgrp/phenotypes/lifespan/female_28c_mean_pheno.txt", header=FALSE, sep="\t")
dat_M_18 <- read.table("/data2/morgante_lab/data/dgrp/phenotypes/lifespan/male_18c_mean_pheno.txt", header=FALSE, sep="\t")
dat_M_25 <- read.table("/data2/morgante_lab/data/dgrp/phenotypes/lifespan/male_25c_mean_pheno.txt", header=FALSE, sep="\t")
dat_M_28 <- read.table("/data2/morgante_lab/data/dgrp/phenotypes/lifespan/male_28c_mean_pheno.txt", header=FALSE, sep="\t")
colnames(dat_F_18) <- colnames(dat_F_25) <- colnames(dat_F_28) <- colnames(dat_M_18) <- colnames(dat_M_25) <- colnames(dat_M_28) <- 
  c("line_id", "y")
rownames(dat_F_18) <- dat_F_18$line_id
rownames(dat_F_25) <- dat_F_25$line_id
rownames(dat_F_28) <- dat_F_28$line_id
rownames(dat_M_18) <- dat_M_18$line_id
rownames(dat_M_25) <- dat_M_25$line_id
rownames(dat_M_28) <- dat_M_28$line_id

###Remove missing phenotypes
dat_F_18 <- dat_F_18[complete.cases(dat_F_18), ]
dat_F_25 <- dat_F_25[complete.cases(dat_F_25), ]
dat_F_28 <- dat_F_28[complete.cases(dat_F_28), ]
dat_M_18 <- dat_M_18[complete.cases(dat_M_18), ]
dat_M_18 <- dat_M_18[complete.cases(dat_M_18), ]
dat_M_18 <- dat_M_18[complete.cases(dat_M_18), ]

###Add environmeantal variables
dat_F_18$temp <- 18
dat_F_18$sex <- 1
dat_F_25$temp <- 25
dat_F_25$sex <- 1
dat_F_28$temp <- 28
dat_F_28$sex <- 1
dat_M_18$temp <- 18
dat_M_18$sex <- 0
dat_M_25$temp <- 25
dat_M_25$sex <- 0
dat_M_28$temp <- 28
dat_M_28$sex <- 0

###Get IDs of lines without missing values
id_count <- table(c(dat_F_18$line_id, dat_F_25$line_id, dat_F_28$line_id, dat_M_18$line_id, dat_M_25$line_id, dat_M_28$line_id))
ids_complete_cases <- names(which(id_count == 6))

###Merge datasets
dat <- rbind(dat_F_18[ids_complete_cases, ],
             dat_F_25[ids_complete_cases, ],
             dat_F_28[ids_complete_cases, ],
             dat_M_18[ids_complete_cases, ],
             dat_M_25[ids_complete_cases, ],
             dat_M_28[ids_complete_cases, ])
rownames(dat) <- 1:nrow(dat)

###Create id of the observation
dat$obs_id <- paste0(dat$line_id, "_", dat$temp, "_", dat$sex)

###Save data
saveRDS(dat, file=output_data)

###Write out list of lines for plink
write.table(cbind(ids_complete_cases, ids_complete_cases), file=output_lines_ids, col.names = FALSE,
            row.names = FALSE, quote=FALSE, sep="\t")
