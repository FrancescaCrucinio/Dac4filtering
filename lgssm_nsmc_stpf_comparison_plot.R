# read data
d <- 32
df <- read.csv("data/lgssm/lgssm_d32N100ID1")
df$N <- "10^2"
for (id in 2:50){
  filename <- paste("data/lgssm/lgssm_d32N100ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^2"
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  filename <- paste("data/lgssm/lgssm_d32N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^3"
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  filename <- paste("data/lgssm/lgssm_d32N10000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^4"
  df <- rbind(df, dfnew)
}
df <- df[, -1]
Time.step <- ncol(df) - 5
colnames(df)[(Time.step+1):(Time.step+3)] <- c("w1", "ks", "runtime")

# distances
distances <- data.frame(aggregate(w1 ~ algo + runtime + N, data = df, FUN = "mean"))
distances$ks <- aggregate(ks ~ algo + runtime + N, data = df, FUN = "mean")$ks
time_means <- aggregate(runtime ~ algo + N, data = df, FUN= "mean" )
time_means <- time_means[order(time_means$algo, time_means$N), ]
distances <- distances[order(distances$algo, distances$N), ]
distances$runtime_mean <- rep(time_means$runtime, each = 1)
# Wasserstein-1
ggplot(data = distances, aes(x = runtime_mean, y = w1, group = interaction(algo, N), fill = algo, colour = algo)) +
  geom_boxplot(coef = 2) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15))
# ggsave("lgssm32_w1.pdf", width = 10, height = 8, dpi = 300)
# Kolmogorov-Smirnov
ggplot(data = distances, aes(x = runtime_mean, y = ks, group = interaction(algo, N), fill = algo, colour = algo)) +
  geom_boxplot(coef = 2) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15))
# ggsave("lgssm32_ks.pdf", width = 10, height = 8, dpi = 300)
# RMSE
tmp <- aggregate(. ~ algo + N, data = df, FUN = "mean")
rmse_data <- data.frame(rep(1:Time.step, times = 3), rep(tmp$algo, each = 100), rep(tmp$N, each = 100))
colnames(rmse_data) <- c("Time.step", "algo", "N")
rmse_data <- rmse_data[order(rmse_data$algo, rmse_data$N), ]
tmp <- tmp[order(tmp$algo, tmp$N), ]
rmse_data$rmse <- as.vector(t(as.matrix(tmp[, 3:102])))
ggplot(data = rmse_data, aes(x = Time.step, y = rmse, group = algo, colour = algo)) +
  geom_line(size = 2) +
  scale_y_continuous(trans='log10') +
  facet_wrap(~N, ncol = 2, nrow = 2) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15))
# ggsave("lgssm32_rmse.pdf", width = 12, height = 8, dpi = 300)
