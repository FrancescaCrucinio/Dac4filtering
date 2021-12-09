# read data
df <- read.csv("data/resampling_comparison_d8N1000ID1")
for (id in 2:50){
  filename <- paste("data/resampling_comparison_d8N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  filename <- paste("data/resampling_comparison_d8N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  filename <- paste("data/resampling_comparison_d8N10000ID", id, sep = "")
  dfnew <- read.csv(filename)
  df <- rbind(df, dfnew)
}
df <- df[, -1]
d <- ncol(df) - 3
df$d_means <- rowMeans(df[, 1:d])
df$N <- as.factor(rep(c("10^3"), each = 50*8))
df <- df[order(df$algo),]
time_means <- aggregate(elapsed ~ algo + N + mutation, data = df, FUN= "mean" )
time_means <- time_means[order(time_means$algo), ]
df$runtime <- rep(time_means$elapsed, each = 50)
# time
ggplot(data = df, aes(x = runtime, y = d_means, group = interaction(algo, mutation), fill = algo, colour = algo)) +
  geom_boxplot(aes(x = runtime, y = d_means, alpha = as.factor(mutation)), coef = 6) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  guides(alpha = FALSE) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=15))

df <- read.csv("data/adaptive_lgssm.csv", col.names = c("u", "m"))
df$u <- as.factor(df$u)
# histogram of m
ggplot(data = df, aes(x = m)) +
  geom_histogram() +
  facet_grid(~u) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=15))
ggsave("histogram_m_d8_N1000.pdf")
