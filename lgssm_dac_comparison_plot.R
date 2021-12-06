# read data
df <- read.csv("data/resampling_comparison_d8N100ID1")
for (id in 2:50){
  filename <- paste("data/resampling_comparison_d8N100ID", id, sep = "")
  dfnew <- read.csv(filename)
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  filename <- paste("data/resampling_comparison_d8N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  df <- rbind(df, dfnew)
}
df <- df[, -1]
d <- ncol(df) - 2
df$d_means <- rowMeans(df[, 1:d])
df$N <- as.factor(rep(c("10^2", "10^3"), each = 50*8))
df <- df[order(df$algo),]
time_means <- aggregate(elapsed ~ algo + N, data = df, FUN= "mean" )
time_means <- time_means[order(time_means$algo), ]
df$runtime <- rep(time_means$elapsed, each = 50)
# time
ggplot(data = df, aes(x = runtime, y = d_means, group = interaction(algo, N), fill = algo, N)) +
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
