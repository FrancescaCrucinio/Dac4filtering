# read data
df <- read.csv("data/resampling_comparison_d8N1000ID2")
for (id in 2:10){
  filename <- paste("data/resampling_comparison_d8N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  df <- rbind(df, dfnew)
}
for (id in 1:10){
  filename <- paste("data/resampling_comparison_d8N10000ID", id, sep = "")
  dfnew <- read.csv(filename)
  df <- rbind(df, dfnew)
}
d <- ncol(df) - 4
df$d_means <- rowMeans(df[, 2:(d+1)])
# time
ggplot(data = df, aes(x = runtime, y = d_means, group = algo, fill = algo)) +
  geom_boxplot(aes(x = runtime, y = d_means), coef = 6, width = 10) +
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
