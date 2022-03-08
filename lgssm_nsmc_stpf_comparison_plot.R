library(ggpubr)
# read data
d <- 256
df <- read.csv(paste0("data/lgssm_tempering/lgssm_d", d, "N100ID1"))
df$N <- "10^2"
df$run <- 1
for (id in 2:50){
  filename <- paste0("data/lgssm_tempering/lgssm_d", d, "N100ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^2"
  dfnew$run <- id
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  filename <- paste0("data/lgssm_tempering/lgssm_d", d, "N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^3"
  dfnew$run <- id
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  filename <- paste("data//lgssm_tempering/lgssm_d32N10000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^4"
  dfnew$run <- id
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  filename <- paste("data/lgssm/v2_lgssm_d32N1e+05ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^5"
  dfnew$run <- id
  df <- rbind(df, dfnew)
}
df <- df[, -1]
df <- df[df$algo != "dac_ada", ]
Time.step <- ncol(df) - 6
colnames(df)[(Time.step+1):(Time.step+3)] <- c("w1", "ks", "runtime")

# distances
distances <- data.frame(aggregate(w1 ~ algo + runtime + N + run, data = df, FUN = "mean"))
distances$ks <- aggregate(ks ~ algo + runtime + N + run, data = df, FUN = "mean")$ks
time_means <- aggregate(runtime ~ algo + N, data = df, FUN= "mean" )
time_means <- time_means[order(time_means$algo, time_means$N), ]
distances <- distances[order(distances$algo, distances$N), ]
distances$runtime_mean <- rep(time_means$runtime, each = 50)
distances_mean <- data.frame(aggregate(w1 ~ algo + N, data = distances, FUN = "mean"))
distances_mean$ks <- aggregate(ks ~ algo + N, data = distances, FUN = "mean")$ks
distances_mean <- merge(distances_mean, time_means, by=c("algo", "N"))
# Wasserstein-1
# w1_plot <-
ggplot(data = distances, aes(x = runtime_mean, y = w1, group = interaction(algo, N), fill = algo, colour = algo)) +
  geom_boxplot(coef = 6, width = 0.1, alpha = 0.1, lwd = 1) +
  geom_point(data = distances_mean, shape = 4, lwd = 1, aes(x = runtime, y = w1, group = interaction(algo, N), fill = algo, colour = algo)) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15), legend.position="none")
# my_legend <- get_legend(w1_plot)
# as_ggplot(my_legend)
# ggsave("lgssm256_w1.pdf", width = 10, height = 8, dpi = 300)
# ggsave("lgssm32_legend.pdf", width = 10, height = 8, dpi = 300)
# Kolmogorov-Smirnov
ggplot(data = distances, aes(x = runtime_mean, y = ks, group = interaction(algo, N), fill = algo, colour = algo)) +
  geom_boxplot(coef = 10, width = 0.1, alpha = 0.1, lwd = 1) +
  geom_point(data = distances_mean, shape = 4, lwd = 1, aes(x = runtime, y = ks, group = interaction(algo, N), fill = algo, colour = algo)) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=25),
        text = element_text(size=20), legend.position="none")
# ggsave("lgssm256_ks.pdf", width = 10, height = 8, dpi = 300)
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
# ggsave("lgssm256_rmse.pdf", width = 12, height = 8, dpi = 300)
