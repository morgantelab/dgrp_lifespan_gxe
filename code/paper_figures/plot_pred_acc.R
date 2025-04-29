###Load libraries
library(ggplot2)
library(cowplot)

set.seed(1)

reps <- c(1,2,4,5,7,9)

###Random lines
model <- c('G_het_var','E_het_var','GandE_het_var','GxE_het_var','mvgblup','rrm')

i <- 0

n_col <- 3
n_row <- length(reps) * length(model)
res <- as.data.frame(matrix(NA, ncol=n_col, nrow=n_row))
colnames(res) <- c("rep", "model", "r2")

for(met in model){
  for(repp in reps){
    dat <- readRDS(paste0("../output/", met, "_fit/dgrp_lifespan_gxe_original_pheno_random_lines_", met, "_fit_", repp, ".rds"))
      
    i <- i + 1
      
    res[i, 1] <- repp
    res[i, 2] <- met
    if(met == "mvgblup"){
      y_test <- do.call("c", dat$y_test)
      yhat_test <- do.call("c", dat$yhat_test)
      res[i, 3] <- summary(lm(y_test ~ yhat_test))$r.squared
    } else {
      res[i, 3] <- summary(lm(dat$y_test ~ dat$yhat_test))$r.squared
    }
  }
}

res <- transform(res,
                 model=factor(model, levels=c("G_het_var", "E_het_var", "GandE_het_var", "GxE_het_var", "mvgblup", "rrm"),
                              labels=c("G-BLUP", "E-BLUP", "GE-BLUP", "GxE-BLUP", "mvG-BLUP", "RRM")))

p_rl <- ggplot(res, aes(x = model, y = r2, fill = model)) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  stat_summary(fun=mean, geom="point", shape=23,
               position = position_dodge2(width = 0.87,   
                                          preserve = "single")) +
  scale_fill_manual(values = c("pink", "red", "yellow", "orange", "green", "blue")) +
  ylim(0, 0.91) +
  labs(x = "Model", y = expression(italic(R)^2), fill="Method", title="") +
  theme_cowplot(font_size = 18) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


###Random observations
model <- c('G_het_var','E_het_var','GandE_het_var','GxE_het_var','mvgblup','rrm')

i <- 0

n_col <- 3
n_row <- length(reps) * length(model)
res <- as.data.frame(matrix(NA, ncol=n_col, nrow=n_row))
colnames(res) <- c("rep", "model", "r2")

for(met in model){
  for(repp in reps){
    dat <- readRDS(paste0("../output/", met, "_fit/dgrp_lifespan_gxe_original_pheno_random_obs_", met, "_fit_", repp, ".rds"))
    
    i <- i + 1
    
    res[i, 1] <- repp
    res[i, 2] <- met
    if(met == "mvgblup"){
      y_test <- do.call("c", dat$y_test)
      yhat_test <- do.call("c", dat$yhat_test)
      res[i, 3] <- summary(lm(y_test ~ yhat_test))$r.squared
    } else {
      res[i, 3] <- summary(lm(dat$y_test ~ dat$yhat_test))$r.squared
    }
  }
}

res <- transform(res,
                 model=factor(model, levels=c("G_het_var", "E_het_var", "GandE_het_var", "GxE_het_var", "mvgblup", "rrm"),
                              labels=c("G-BLUP", "E-BLUP", "GE-BLUP", "GxE-BLUP", "mvG-BLUP", "RRM")))

p_ro <- ggplot(res, aes(x = model, y = r2, fill = model)) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  stat_summary(fun=mean, geom="point", shape=23,
               position = position_dodge2(width = 0.87,   
                                          preserve = "single")) +
  scale_fill_manual(values = c("pink", "red", "yellow", "orange", "green", "blue")) +
  ylim(0, 0.91) +
  labs(x = "Model", y = expression(italic(R)^2), fill="Method", title="") +
  theme_cowplot(font_size = 18) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


###New environment
n_reps <- 6
model <- c('G_het_var','E_het_var','GandE_het_var','GxE_het_var','rrm')

i <- 0

n_col <- 3
n_row <- length(reps) * length(model)
res <- as.data.frame(matrix(NA, ncol=n_col, nrow=n_row))
colnames(res) <- c("rep", "model", "r2")

for(met in model){
  for(repp in 1:n_reps){
    dat <- readRDS(paste0("../output/", met, "_fit/dgrp_lifespan_gxe_original_pheno_temp_sex_", met, "_fit_", repp, ".rds"))
    
    i <- i + 1
    
    res[i, 1] <- repp
    res[i, 2] <- met
    if(met == "mvgblup"){
      y_test <- do.call("c", dat$y_test)
      yhat_test <- do.call("c", dat$yhat_test)
      res[i, 3] <- summary(lm(y_test ~ yhat_test))$r.squared
    } else {
      res[i, 3] <- summary(lm(dat$y_test ~ dat$yhat_test))$r.squared
    }
  }
}

res <- transform(res,
                 model=factor(model, levels=c("G_het_var", "E_het_var", "GandE_het_var", "GxE_het_var", "rrm"),
                              labels=c("G-BLUP", "E-BLUP", "GE-BLUP", "GxE-BLUP", "RRM")))

p_ne <- ggplot(res, aes(x = model, y = r2, fill = model)) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  stat_summary(fun=mean, geom="point", shape=23,
               position = position_dodge2(width = 0.87,   
                                          preserve = "single")) +
  scale_fill_manual(values = c("pink", "red", "yellow", "orange", "blue")) +
  ylim(0, 0.91) +
  labs(x = "Model", y = expression(italic(R)^2), fill="Method", title="") +
  theme_cowplot(font_size = 18) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


###Combine plots
p <- plot_grid(p_rl,
               p_ro,
               p_ne,
               nrow=1,
               labels = c("A", "B", "C"))

ggsave("../analysis/paper_figures/Fig3.eps", plot=p, device="eps", units="in", height=5, width=10)
