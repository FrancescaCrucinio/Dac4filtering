library(plyr)
# linear Gaussian
df <- read.csv("data/adaptive_resampling/adaptive_lgssm_d32N1000T100.csv", col.names = c("u", "m"))
ggplot(data = df, aes(x = m)) +
  geom_histogram() +
  facet_grid(~factor(u, levels=c('1','2','3','4', '5'))) +
  scale_y_continuous(trans='log1p') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Linear Gaussian")
# ggsave("m_adaptive32_lgssm.pdf", width = 10, height = 5, dpi = 300)

# Spatial
df <- read.csv("data/adaptive_resampling/adaptive_nl.csv", col.names = c("u", "direction", "m"))
# histogram of m
ggplot(data = df, aes(x = m)) +
  geom_histogram() +
  facet_grid(~interaction(factor(u, levels = c('1', '2', '3')), factor(direction, levels = c('h', 'v')))) +
  scale_y_continuous(trans='log1p') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Non-linear Spatial")
# ggsave("m_adaptive8x8_nl.pdf", width = 10, height = 5, dpi = 300)

df <- read.csv("data/adaptive_resampling/adaptive_nl_cov.csv", col.names = c("u", "direction", "m"))
# histogram of m
ggplot(data = df, aes(x = m)) +
  geom_histogram() +
  facet_grid(~interaction(factor(u, levels = c('1', '2', '3')), factor(direction, levels = c('h', 'v')))) +
  scale_y_continuous(trans='log1p') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Non-linear Spatial w/ covariance")
# ggsave("m_adaptive8x8_nl_cov.pdf", width = 10, height = 5, dpi = 300)
