library(plyr)
# df <- read.csv("data/adaptive_resampling/adaptive_nl_d64N1000T10.csv", col.names = c("u", "direction", "m"))
df <- read.csv("data/adaptive_resampling/adaptive_lgssm_d32N1000T100.csv", col.names = c("u", "m"))
df$u <- as.factor(mapvalues(df$u, from=c(1, 2, 3, 4, 5), to=c(5, 4, 3, 2, 1)))
# histogram of m
ggplot(data = df, aes(x = m)) +
  geom_histogram() +
  facet_grid(~interaction(factor(u, levels=c('5','4','3','2', '1')))) +
  scale_y_continuous(trans='log1p') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Linear Gaussian")
# ggsave("m_adaptive32_lgssm.pdf", width = 10, height = 5, dpi = 300)
