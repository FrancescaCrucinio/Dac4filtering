# Plot of W1 and KS distance for linear Gaussian SSM
library(ggpubr)
# read data
df <- read.csv(paste0("data/lgssm_tempering/lgssm_d32N100ID1"))
df$N <- "10^2"
df$d <- "32"
df$run <- 1
dfnew <- read.csv(paste0("data/lgssm_tempering/lgssm_d256N100ID1"))
dfnew$N <- "10^2"
dfnew$d <- "256"
dfnew$run <- 1
df <- rbind(df, dfnew)
dfnew <- read.csv(paste0("data/lgssm_tempering/lgssm_d2048N100ID1"))
dfnew$N <- "10^2"
dfnew$d <- "2048"
dfnew$run <- 1
df <- rbind(df, dfnew)
for (id in 2:50){
  filename <- paste0("data/lgssm_tempering/lgssm_d32N100ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^2"
  dfnew$d <- "32"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/lgssm_d256N100ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^2"
  dfnew$d <- "256"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/lgssm_d2048N100ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^2"
  dfnew$d <- "2048"
  dfnew$run <- id
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  filename <- paste0("data/lgssm_tempering/lgssm_d32N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^3"
  dfnew$d <- "32"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/lgssm_d256N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^3"
  dfnew$d <- "256"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/lgssm_d2048N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^3"
  dfnew$d <- "2048"
  dfnew$run <- id
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  filename <- paste0("data/lgssm_tempering/lgssm_d32N10000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^4"
  dfnew$d <- "32"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/lgssm_d256N10000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^4"
  dfnew$d <- "256"
  dfnew$run <- id
  df <- rbind(df, dfnew)
}
df <- df[, -1]
df <- df[df$algo != "dac_ada", ]
Time.step <- ncol(df) - 7
colnames(df)[(Time.step+1):(Time.step+3)] <- c("w1", "ks", "runtime")


# distances
distances <- data.frame(aggregate(w1 ~ algo + runtime + N + run + d, data = df, FUN = "mean"))
distances$ks <- aggregate(ks ~ algo + runtime + N + run + d, data = df, FUN = "mean")$ks
time_means <- aggregate(runtime ~ algo + N + d, data = df, FUN= "mean" )
time_means <- time_means[order(time_means$algo, time_means$N, time_means$d), ]
distances <- distances[order(distances$algo, distances$N, distances$d), ]
distances$runtime_mean <- rep(time_means$runtime, each = 50)
distances_mean <- data.frame(aggregate(w1 ~ algo + N + d, data = distances, FUN = "mean"))
distances_mean$ks <- aggregate(ks ~ algo + N + d, data = distances, FUN = "mean")$ks
distances_mean <- merge(distances_mean, time_means, by=c("algo", "N", "d"))
# Wasserstein-1
ggplot(data = distances, aes(x = runtime_mean, y = w1, group = interaction(algo, N, d), fill = algo, colour = algo)) +
  geom_boxplot(coef = 6, width = 0.1, alpha = 0.1, lwd = 1) +
  geom_point(data = distances_mean, shape = 4, lwd = 1, aes(x = runtime, y = w1, group = interaction(algo, N, d), fill = algo, colour = algo)) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~factor(d, levels = c("32", "256", "2048")), ncol = 1, nrow = 3) +
  scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF")) +
  scale_color_manual(values=c("#F8766D", "#00BA38", "#619CFF")) +
  theme(axis.title.x=element_blank(),  axis.text = element_text(size=30),
        axis.title.y=element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30), legend.position="none", panel.spacing = unit(2, "lines"))
# ggsave("lgssm_w1.pdf", width = 20, height = 30, dpi = 300)
# Kolmogorov-Smirnov
ggplot(data = distances, aes(x = runtime_mean, y = ks, group = interaction(algo, N, d), fill = algo, colour = algo)) +
  geom_boxplot(coef = 10, width = 0.1, alpha = 0.1, lwd = 1) +
  geom_point(data = distances_mean, shape = 4, lwd = 1, aes(x = runtime, y = ks, group = interaction(algo, N, d), fill = algo, colour = algo)) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~factor(d, levels = c("32", "256", "2048")), ncol = 1, nrow = 3) +
  scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF")) +
  scale_color_manual(values=c("#F8766D", "#00BA38", "#619CFF")) +
  theme(axis.title.x=element_blank(), axis.text = element_text(size=30),
        axis.title.y=element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30), legend.position="none", panel.spacing = unit(2, "lines"))
# ggsave("lgssm_ks.pdf", width = 20, height = 30, dpi = 300)
