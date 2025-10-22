###Load libraries needed
library(optparse)
library(tidyr)
library(dplyr)
library(sommer)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--grm_file"), type="character")
parser <- add_option(parser, c("--pheno_file"), type="character")
parser <- add_option(parser, c("--cv_file"), type="character")
parser <- add_option(parser, c("--rep"), type="integer")
parser <- add_option(parser, c("--fixed_str"), type="character")
parser <- add_option(parser, c("--random_str"), type="character")
parser <- add_option(parser, c("--residual_str"), type="character")
parser <- add_option(parser, c("--verbose"), type="logical")
parser <- add_option(parser, c("--output_file"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

grm_file <- outparse$grm_file
pheno_file <- outparse$pheno_file
cv_file <- outparse$cv_file
repp <- outparse$rep
fixed_str <- outparse$fixed_str
random_str <- outparse$random_str
residual_str <- outparse$residual_str
verbose <- outparse$verbose
output_file <- outparse$output_file
seed <- outparse$seed

###Set seed
set.seed(seed)

###Load data
G <- readRDS(grm_file)
dat <- readRDS(pheno_file)
cv <- readRDS(cv_file)

###Compute environment means
env_means <- dat %>% group_by(sex, temp) %>% summarise(avg_y_by_env = mean(y)) %>% as.data.frame()
env_means <- env_means[order(env_means$avg_y_by_env),]

###Add environment means to main data
dat <- left_join(dat, env_means, by=join_by(sex,temp))

###Order GRM according to phenotype data
G <- G[unique(dat$line_id), unique(dat$line_id)]

###Set test set observations to missing
dat$yNa <- dat$y
whichNa <- which(dat$obs_id %in% cv[[repp]])
dat$yNa[whichNa] <- NA

###Fit random regression model
fit <- mmer(fixed=as.formula(fixed_str),
             random=as.formula(random_str), 
             rcov=as.formula(residual_str),
             nIters=100,
             method="NR",
             dateWarning=FALSE, 
             date.warning=FALSE, 
             verbose=verbose,
             data=dat)

###Extract random regression GEBVs
leg_max_order <- as.numeric(substring(gsub("\\s", "", unlist(strsplit(random_str, ","))[2]), 1, 1))
RnReg <- vector("list", leg_max_order)

for(i in 0:leg_max_order){
  RnReg[[i + 1]] <- unlist(eval(parse(text=paste0("fit$U$`leg", i, ":line_id`"))))
}

RnReg <- do.call("rbind", RnReg)
colnames(RnReg) <- unique(dat$line_id)
envv <- unique(dat$avg_y_by_env)
Phi <- leg(envv, leg_max_order)
dat$env <- paste(dat$sex, dat$temp, sep="_")
rownames(Phi) <- unique(dat$env)
GEBV <- t(Phi %*% RnReg)

###Compute predictions
if(grepl("leg", fixed_str)){
  leg_max_order_fix <- as.numeric(substring(gsub("\\s", "", unlist(strsplit(fixed_str, ","))[2]), 1, 1))
  Phi_f <- leg(envv, leg_max_order_fix)
  mean_function <- drop(Phi_f %*% fit$Beta$Estimate)
  yHat <- t(mean_function + t(GEBV))
} else {
  effects_sel <- gsub("\\s", "", unlist(strsplit(unlist(strsplit(fixed_str, "~"))[2], "+", fixed=TRUE)))
  
  X <- model.matrix(as.formula(paste("y ~", paste(effects_sel, collapse = "+"))), data=dat)
  dat$Xbeta <- X %*% fit$Beta$Estimate

  mean_function <- dat %>%
    dplyr::select(line_id, all_of(effects_sel), Xbeta) %>%
    pivot_wider(names_from = all_of(effects_sel), values_from = Xbeta) %>% 
    dplyr::select(!line_id) %>%
    as.matrix()
  
  rownames(mean_function) <- unique(dat$line_id)
  
  yHat <- mean_function + GEBV
}

###Extract predictions for the test set
yHat_test <- vector("numeric", length(whichNa)) 

it <- 0

for(i in 1:nrow(dat)){
  if(is.na(dat$yNa[i])){
    it <- it+1
    row_to_keep <- which(rownames(yHat) == dat$line_id[i]) 
    col_to_keep <- which(colnames(yHat) == dat$env[i])
    yHat_test[it] <- yHat[row_to_keep, col_to_keep]
    names(yHat_test)[it] <- dat$obs_id[i]
  }
}

###Compile results and write them to a file
res <- list(model_fit=fit, y_test=dat$y[whichNa], yhat_test=yHat_test)
saveRDS(res, output_file)



