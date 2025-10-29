###Load libraries needed
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

dat <- readRDS("/data2/morgante_lab/fabiom/dgrp_lifespan_gxe/data/dgrp_lifespan_gxe_original_pheno_env.rds")

###Compute environment means
env_means <- dat %>% group_by(sex, temp) %>% summarise(avg_y_by_env = mean(y)) %>% as.data.frame()
env_means <- env_means[order(env_means$avg_y_by_env),]

###Add environment means to main data
dat <- left_join(dat, env_means, by=join_by(sex,temp))

##Generate a palette
my_brewer_palette_function <- colorRampPalette(brewer.pal(11, "Set3"))
my_palette <- my_brewer_palette_function(186)

p <- ggplot(dat, aes(x=avg_y_by_env, y=y, group=line_id, color=line_id)) +
  geom_line() +
  scale_color_manual(values = my_palette) +
  labs(x = "Environmental value", y = "Phenotypic value", title="") +
  theme_cowplot(font_size = 18) +
  theme(legend.position="none")

ggsave("../analysis/paper_figures/FigS1.eps", plot=p, device="eps", units="in", height=5, width=7)


