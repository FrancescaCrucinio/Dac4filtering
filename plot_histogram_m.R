library(plyr)
library(ggplot2)
# linear Gaussian
df <- read.csv("data/adaptive_resampling/adaptive_lgssm.csv", col.names = c("u", "m"))
ggplot(data = df, aes(x = m)) +
  geom_histogram() +
  facet_grid(~factor(u, levels=c('1','2','3','4', '5'))) +
  scale_y_continuous(trans='log1p') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Linear Gaussian")
# ggsave("marginal_adaptive32_lgssm.pdf", width = 10, height = 5, dpi = 300)

# Spatial
df <- read.csv("data/adaptive_resampling/adaptive_spatial.csv", col.names = c("u", "direction", "m"))
# histogram of m
ggplot(data = df, aes(x = m)) +
  geom_histogram() +
  facet_wrap(~interaction(factor(u, levels = c('1', '2', '3', '4', '5')), factor(direction, levels = c('h', 'v'))), nrow = 2, ncol = 5) +
  scale_y_continuous(trans='log1p') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=25), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Spatial")
# ggsave("m_adaptive_spatial.pdf", width = 14, height = 12, dpi = 300)
