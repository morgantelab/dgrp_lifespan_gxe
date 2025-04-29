###Set seed
set.seed(1)

###Load data
dat <- readRDS("../data/dgrp_lifespan_gxe_original_pheno_env.rds")

###Extract environmental variables and scale them
X <- dat[, c("temp", "sex")]
We <- scale(X)

###Fit linear regression model
##Scaled X
fit <- lm(dat$y ~ We)
summary(fit)

##Original X
fit1 <- lm(y ~ temp + sex, data=dat)
summary(fit1)
