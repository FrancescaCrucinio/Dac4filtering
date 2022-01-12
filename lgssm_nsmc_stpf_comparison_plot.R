# read data
df <- read.csv("data/lgssm_d32N100ID1")
for (id in 2:50){
  filename <- paste("data/lgssm_d32N100ID", id, sep = "")
  dfnew <- read.csv(filename)
  df <- rbind(df, dfnew)
}
df <- df[, -1]
d <- nrow(df)/(3*2)
Time.step <- ncol(df) - 4
colnames(df)[(Time.step+1):(Time.step+3)] <- c("w1", "ks", "runtime")

# distances
distances <- data.frame(aggregate(w1 ~ algo + runtime, data = df, FUN = "mean"))
distances$ks <- aggregate(ks ~ algo + runtime, data = df, FUN = "mean")$ks
time_means <- aggregate(runtime ~ algo , data = df, FUN= "mean" )
time_means <- time_means[order(time_means$algo), ]
distances <- distances[order(distances$algo), ]
distances$runtime_mean <- rep(time_means$runtime, each = 2)
# Wasserstein-1
ggplot(data = distances, aes(x = runtime_mean, y = w1, group = interaction(algo), fill = algo, colour = algo)) +
  geom_boxplot(aes(x = runtime, y = w1), coef = 2) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15))
# Kolmogorov-Smirnov
ggplot(data = distances, aes(x = runtime_mean, y = ks, group = interaction(algo), fill = algo, colour = algo)) +
  geom_boxplot(aes(x = runtime, y = ks), coef = 2) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15))
# RMSE
tmp <- aggregate(. ~ algo, data = df, FUN = "mean")
rmse_data <- data.frame(rep(1:Time.step, times = 3), rep(tmp$algo, each = 100))
colnames(rmse_data) <- c("Time.step", "algo")
rmse_data$rmse <- as.vector(t(as.matrix(tmp[, 2:101])))
ggplot(data = rmse_data, aes(x = Time.step, y = rmse, group = algo, colour = algo)) +
  geom_line(aes(x = Time.step, y = rmse)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15))
