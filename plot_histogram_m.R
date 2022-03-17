library(plyr)
df <- read.csv("data/adaptive_nl.csv", col.names = c("u", "direction", "m"))
# df$u <- as.factor(mapvalues(df$u, from=c(1, 2, 3, 4, 5), to=c(5, 4, 3, 2, 1)))
# histogram of m
ggplot(data = df, aes(x = m)) +
  geom_histogram() +
  facet_grid(~interaction(factor(u, levels=c('5','4','3','2', '1')), direction)) +
  scale_y_continuous(trans='log1p') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=15),
        text = element_text(size=15))
# ggsave("m_adaptive32_car.pdf", width = 10, height = 5, dpi = 300)
