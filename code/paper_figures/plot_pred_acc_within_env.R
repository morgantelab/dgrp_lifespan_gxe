###Load libraries
library(ggplot2)
library(cowplot)
library(dplyr)

set.seed(1)

reps <- c(1,2,4,5,7,9)

###Random lines
model <- c('G','E','GandE','GxE','mvgblup','rrm')

i <- 0

n_col <- 4
n_row <- length(reps) * length(model) * 6
res <- as.data.frame(matrix(NA, ncol=n_col, nrow=n_row))
colnames(res) <- c("rep", "env", "model", "r2")

for(met in model){
  for(repp in reps){
    if(met %in% c('G','E','GandE','GxE')){
      dat <- readRDS(paste0("../output/", met, "_fit/dgrp_lifespan_gxe_original_pheno_random_lines_", met, "_reml_fit_", repp, ".rds"))
    } else {
      dat <- readRDS(paste0("../output/", met, "_fit/dgrp_lifespan_gxe_original_pheno_random_lines_", met, "_fit_", repp, ".rds"))
      
      if(met=="mvgblup"){
        for(t in 1:6){
          temp <- unlist(strsplit(colnames(dat$model_fit$ETAHat)[t], "_"))[2]
          sex <- unlist(strsplit(colnames(dat$model_fit$ETAHat)[t], "_"))[3]
          names(dat$y_test[[t]]) <- paste(names(dat$y_test[[t]]), temp, sex, sep="_")
          names(dat$yhat_test[[t]]) <- paste(names(dat$yhat_test[[t]]), temp, sex, sep="_")
        }
        
        dat$y_test <- do.call("c", dat$y_test)
        dat$yhat_test <- do.call("c", dat$yhat_test)
      }
      
    }
    
    dat_pheno <- data.frame(y=dat$y_test, yhat=dat$yhat_test, obs=names(dat$yhat_test), temp=as.numeric(NA), sex=as.numeric(NA))
    
    for(k in 1:nrow(dat_pheno)){
      obs <- dat_pheno[k, "obs"]
      obs_splt <- unlist(strsplit(obs, "_"))
      dat_pheno[k, "temp"] <- obs_splt[3]
      dat_pheno[k, "sex"] <- obs_splt[4]
    }
    
    if(i == 0){
      i <- i + 1
    } else {
      i <- i + 6
    }
    
    
    models <- dat_pheno %>%
      group_by(temp, sex) %>%
      do(r2 = summary(lm(y ~ yhat, data = .))$r.squared) %>%
      as.data.frame()
      
    models$env <- paste(models$temp, models$sex, sep="_") 
      
    for(j in 0:5){
      res[i+j, 1] <- repp
      res[i+j, 2] <- models[j+1, "env"]
      res[i+j, 3] <- met
      res[i+j, 4] <- models[j+1, "r2"]
    }
  }
}

##If want to average across environments
res <- res %>% group_by(rep, model) %>% summarise_at(vars(r2), list(r2 = mean)) %>% as.data.frame()

res <- transform(res,
                 model=factor(model, levels=c("G", "E", "GandE", "GxE", "mvgblup", "rrm"),
                              labels=c("G-BLUP", "E-BLUP", "GE-BLUP", "GxE-BLUP", "mvG-BLUP", "RRM")))

p_rl <- ggplot(res, aes(x = model, y = r2, fill = model)) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  stat_summary(fun=mean, geom="point", shape=23,
               position = position_dodge2(width = 0.87,   
                                          preserve = "single")) +
  scale_fill_manual(values = c("pink", "red", "yellow", "orange", "green", "blue")) +
  ylim(0, 0.11) +
  labs(x = "Model", y = expression(italic(R)^2), fill="Method", title="") +
  theme_cowplot(font_size = 18) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

##If want to separate by environment
# res <- transform(res,
#                  model=factor(model, levels=c("G", "E", "GandE", "GxE", "mvgblup", "rrm"),
#                               labels=c("G-BLUP", "E-BLUP", "GE-BLUP", "GxE-BLUP", "mvG-BLUP", "RRM")),
#                  env=factor(env, levels = c("18_1", "18_0", "25_1", "25_0", "28_1", "28_0"),
#                               labels = c("18_F", "18_M", "25_F", "25_M", "28_F", "28_M")))
#
# p_rl <- ggplot(res, aes(x = env, y = r2, fill = model)) +
#   geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
#   stat_summary(fun=mean, geom="point", shape=23,
#                position = position_dodge2(width = 0.87,   
#                                           preserve = "single")) +
#   scale_fill_manual(values = c("pink", "red", "yellow", "orange", "green", "blue")) +
#   ylim(0, 0.3) +
#   labs(x = "Model", y = expression(italic(R)^2), fill="Method", title="") +
#   theme_cowplot(font_size = 18) +
#   theme(legend.position="none",
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 

###Random observations
model <- c('G','E','GandE','GxE','mvgblup', 'rrm')

i <- 0

n_col <- 4
n_row <- length(reps) * length(model) * 6
res <- as.data.frame(matrix(NA, ncol=n_col, nrow=n_row))
colnames(res) <- c("rep", "env", "model", "r2")

for(met in model){
  for(repp in reps){
    
    if(met %in% c('G','E','GandE','GxE')){
      dat <- readRDS(paste0("../output/", met, "_fit/dgrp_lifespan_gxe_original_pheno_random_obs_", met, "_reml_fit_", repp, ".rds"))
    } else {
      dat <- readRDS(paste0("../output/", met, "_fit/dgrp_lifespan_gxe_original_pheno_random_obs_", met, "_fit_", repp, ".rds"))
      if(met=="mvgblup"){
        for(t in 1:6){
          temp <- unlist(strsplit(colnames(dat$model_fit$ETAHat)[t], "_"))[2]
          sex <- unlist(strsplit(colnames(dat$model_fit$ETAHat)[t], "_"))[3]
          names(dat$y_test[[t]]) <- paste(names(dat$y_test[[t]]), temp, sex, sep="_")
          names(dat$yhat_test[[t]]) <- paste(names(dat$yhat_test[[t]]), temp, sex, sep="_")
        }
        
        dat$y_test <- do.call("c", dat$y_test)
        dat$yhat_test <- do.call("c", dat$yhat_test)
      }
      
    }
    
    dat_pheno <- data.frame(y=dat$y_test, yhat=dat$yhat_test, obs=names(dat$yhat_test), temp=as.numeric(NA), sex=as.numeric(NA))
    
    for(k in 1:nrow(dat_pheno)){
      obs <- dat_pheno[k, "obs"]
      obs_splt <- unlist(strsplit(obs, "_"))
      dat_pheno[k, "temp"] <- obs_splt[3]
      dat_pheno[k, "sex"] <- obs_splt[4]
    }
    
    if(i == 0){
      i <- i + 1
    } else {
      i <- i + 6
    }
    
    
    models <- dat_pheno %>%
      group_by(temp, sex) %>%
      do(r2 = summary(lm(y ~ yhat, data = .))$r.squared) %>%
      as.data.frame()
    
    models$env <- paste(models$temp, models$sex, sep="_") 
    
    for(j in 0:5){
      res[i+j, 1] <- repp
      res[i+j, 2] <- models[j+1, "env"]
      res[i+j, 3] <- met
      res[i+j, 4] <- models[j+1, "r2"]
    }
  }
}

##If want to average across environments
res <- res %>% group_by(rep, model) %>% summarise_at(vars(r2), list(r2 = mean)) %>% as.data.frame()

res <- transform(res,
                 model=factor(model, levels=c("G", "E", "GandE", "GxE", "mvgblup", "rrm"),
                              labels=c("G-BLUP", "E-BLUP", "GE-BLUP", "GxE-BLUP", "mvG-BLUP", "RRM")))

p_ro <- ggplot(res, aes(x = model, y = r2, fill = model)) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  stat_summary(fun=mean, geom="point", shape=23,
               position = position_dodge2(width = 0.87,   
                                          preserve = "single")) +
  scale_fill_manual(values = c("pink", "red", "yellow", "orange", "green", "blue")) +
  ylim(0, 0.65) +
  labs(x = "Model", y = expression(italic(R)^2), fill="Method", title="") +
  theme_cowplot(font_size = 18) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

##If want to separate by environment
# res <- transform(res,
#                  model=factor(model, levels=c("G", "E", "GandE", "GxE", "mvgblup", "rrm"),
#                               labels=c("G-BLUP", "E-BLUP", "GE-BLUP", "GxE-BLUP", "mvG-BLUP", "RRM")),
#                  env=factor(env, levels = c("18_1", "18_0", "25_1", "25_0", "28_1", "28_0"),
#                             labels = c("18_F", "18_M", "25_F", "25_M", "28_F", "28_M")))
# 
# p_ro <- ggplot(res, aes(x = env, y = r2, fill = model)) +
#   geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
#   stat_summary(fun=mean, geom="point", shape=23,
#                position = position_dodge2(width = 0.87,   
#                                           preserve = "single")) +
#   scale_fill_manual(values = c("pink", "red", "yellow", "orange", "green", "blue")) +
#   ylim(0, 0.91) +
#   labs(x = "Model", y = expression(italic(R)^2), fill="Method", title="") +
#   theme_cowplot(font_size = 18) +
#   theme(legend.position="none",
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


###Combine plots
p <- plot_grid(p_rl,
               p_ro,
               nrow=1,
               labels = c("A", "B"))

ggsave("../analysis/paper_figures/FigS2.eps", plot=p, device="eps", units="in", height=5, width=10)



