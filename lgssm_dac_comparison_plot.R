# read data
d <- ncol(df) - 3
df$d_means <- rowMeans(df[, 1:d])
# time
ggplot(data = df, aes(x = runtime, y = d_means, group = algo, fill = algo)) +
  geom_boxplot(aes(x = runtime, y = d_means), coef = 6) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=15))
# memory
ggplot(data = df, aes(x = memory, y = d_means, group = algo, fill = algo)) +
  geom_boxplot(aes(x = memory, y = d_means), coef = 6) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=15))
