###Load libraries
library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)


###Variance components
model <- c("G", "E", "GandE", "GxE")

res_var <- as.data.frame(matrix(NA, ncol=6, nrow=4))
colnames(res_var) <- c("model", "VarG", "VarE", "VarGxE", "VarEps", "VarTot")

i <- 0

for(met in model){
  dat <- readRDS(paste0("../output/", met, "_fit/dgrp_lifespan_gxe_original_pheno_", met, "_fit_whole_data.rds"))
  i <- i + 1
  
  res_var[i, 1] <- met
  if(!is.null(dat$ETA$G$varU)){
    res_var[i, 2] <- dat$ETA$G$varU
  } else {
    res_var[i, 2] <- NA
  }
  
  if(!is.null(dat$ETA$E$varU)){
    res_var[i, 3] <- dat$ETA$E$varU
  } else {
    res_var[i, 3] <- NA
  }
  
  if(!is.null(dat$ETA$GxE$varU)){
    res_var[i, 4] <- dat$ETA$GxE$varU
  } else {
    res_var[i, 4] <- NA
  }
  
  res_var[i, 5] <- dat$varE
  res_var[i, 6] <- sum(res_var[i, 2:5], na.rm=TRUE)
}


###PVE
res_var$propG <- res_var$VarG/res_var$VarTot
res_var$propE <- res_var$VarE/res_var$VarTot
res_var$propGxE <- res_var$VarGxE/res_var$VarTot
res_var$propEps <- res_var$VarEps/res_var$VarTot

res_var_long <- res_var %>% select(model, VarG, VarE, VarGxE, VarEps) %>% 
  gather(value="Var", key="Source", VarG, VarE, VarGxE, VarEps)

res_var_long <- transform(res_var_long, model=factor(model, levels=c("G", "E", "GandE", "GxE"),
                                                     labels=c("G", "E", "G+E", "G+E+GxE")),
                          Source=factor(Source, levels=c("VarEps", "VarGxE", "VarE", "VarG")))



p_varcomps <- ggplot(res_var_long, aes(x = model, y = Var, fill = Source)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("lightgreen", "orange", "pink", "blue"), labels=c(expression(sigma[epsilon]^2), 
                                                                                 expression(sigma[ae]^2), 
                                                                                 expression(sigma[e]^2),
                                                                                 expression(sigma[a]^2))) +
  labs(x = "Model", y = expression(italic(PVE)), fill="Source", title="") +
  theme_cowplot(font_size = 18) +
  theme(legend.position="bottom",
        legend.justification.bottom = "center")

###Genetic correlations
dat <- readRDS("../output/mvgblup_fit/dgrp_lifespan_gxe_original_pheno_mvgblup_fit_whole_data.rds")

corg <- cov2cor(dat$ETA$G$Cov$Omega)
envs <- gsub("y_", "", colnames(dat$ETAHat))
envs <- gsub("_1", "_F", envs)
colnames(corg) <- rownames(corg) <- gsub("_0", "_M", envs)

corg[upper.tri(corg)]<- NA

corg_melt <-reshape2::melt(corg)
colnames(corg_melt) <- c("Env1", "Env2", "Cor")

p_gencor <- ggplot(corg_melt, aes(x = Env1, y = Env2, fill = Cor)) +
  geom_tile(color = "white") +
  geom_text(aes(Env1, Env2, label = round(Cor, 2)), color = "black", size = 4) +
  labs(x = "", y = "", fill=expression(r[A]), title="") +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000",
                       na.value = "white") +
  coord_fixed() +
  theme_cowplot(font_size = 18) +
  theme(legend.position="bottom",
        legend.justification.bottom = "center",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 1))

p <- plot_grid(p_varcomps,
               p_gencor,
               labels = c("A", "B"))

ggsave("../analysis/paper_figures/Fig1.eps", plot=p, device="eps", units="in", height=5, width=10)

