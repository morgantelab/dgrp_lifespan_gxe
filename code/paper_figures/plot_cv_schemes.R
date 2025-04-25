###Load libraries
library(ggplot2)
library(cowplot)
library(reshape2)

set.seed(12)

Y <- matrix(1, nrow=20, ncol=5)
rownames(Y) <- paste0("n_", 20:1)
colnames(Y) <- paste0("e_", 1:5)

###Random lines
Yna <- Y
Yna[c(1,4, 10, 17), ] <- NA

Yna_melt <- melt(Yna)
colnames(Yna_melt) <- c("n", "e", "y")

p_rl <- ggplot(Yna_melt, aes(x = e, y = n, fill = y)) +
  geom_tile(color = "white") +
  labs(x = "", y = "", fill=expression(r[A]), title="") +
  scale_fill_gradient2(low = "#FFFFCC",
                       mid = "#075AFF",
                       high = "#FF0000",
                       na.value = "green") +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  theme_cowplot(font_size = 14) +
  theme(legend.position = "none",
        line = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=0.1))

###Random observations
Yna <- Y

Yna[sample(x=1:100, size=20)] <- NA

Yna_melt <- melt(Yna)
colnames(Yna_melt) <- c("n", "e", "y")

p_ro <- ggplot(Yna_melt, aes(x = e, y = n, fill = y)) +
  geom_tile(color = "white") +
  labs(x = "", y = "", fill=expression(r[A]), title="") +
  scale_fill_gradient2(low = "#FFFFCC",
                       mid = "#075AFF",
                       high = "#FF0000",
                       na.value = "green") +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  theme_cowplot(font_size = 14) +
  theme(legend.position = "none",
        line = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=0.1))

###New environment
Yna <- Y
Yna[, 2] <- NA

Yna_melt <- melt(Yna)
colnames(Yna_melt) <- c("n", "e", "y")

p_ne <- ggplot(Yna_melt, aes(x = e, y = n, fill = y)) +
  geom_tile(color = "white") +
  labs(x = "", y = "", title="") +
  scale_fill_gradient2(low = "#FFFFCC",
                       mid = "#075AFF",
                       high = "#FF0000",
                       na.value = "green") +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  theme_cowplot(font_size = 14) +
  theme(legend.position = "none",
        line = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=0.1))

##Combine the plots
p <- plot_grid(p_rl,
               p_ro,
               p_ne,
               labels = c("A", "B", "C"),
               nrow=1)

ggsave("../analysis/paper_figures/Fig2.eps", plot=p, device="eps", units="in", height=5, width=10)
