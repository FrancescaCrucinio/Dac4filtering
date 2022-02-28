d <- 128
# read data
df <- rbind(read.csv(paste0("data/resampling_tempering/nocrossover_resampling_comparison_d", d, "N100ID1")),
            read.csv(paste0("data/resampling_tempering/crossover_resampling_comparison_d", d, "N100ID1")))
df$N <- "10^2"
for (id in 2:50){
  dfnew <- rbind(read.csv(paste0("data/resampling_tempering/nocrossover_resampling_comparison_d", d, "N100ID", id, sep = "")),
                 read.csv(paste0("data/resampling_tempering/crossover_resampling_comparison_d", d, "N100ID", id, sep = "")))
  dfnew$N <- "10^2"
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  dfnew <- rbind(read.csv(paste0("data/resampling_tempering/nocrossover_resampling_comparison_d", d, "N1000ID", id, sep = "")),
                 read.csv(paste0("data/resampling_tempering/crossover_resampling_comparison_d", d, "N1000ID", id, sep = "")))
  dfnew$N <- "10^3"
  df <- rbind(df, dfnew)
}
df <- df[, -1]
d <- ncol(df) - 4
df$d_means <- rowMeans(df[, 1:d])
df <- df[order(df$algo, df$mutation, df$N),]
time_means <- aggregate(elapsed ~ algo + mutation + N, data = df, FUN = "mean" )
time_means <- time_means[order(time_means$algo, time_means$mutation, time_means$N), ]
df$runtime <- rep(time_means$elapsed, each = 50)
aggregate(d_means ~ algo + mutation + N, data = df, FUN = "mean")
# time
ggplot(data = df, aes(x = runtime, y = d_means, group = interaction(algo, mutation, N), fill = algo, colour = algo)) +
  geom_boxplot(aes(alpha = as.factor(mutation)), coef = 6, width = 20) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  guides(alpha = "none") +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15))
# ggsave("res_time_medium.pdf", width = 10, height = 5, dpi = 300)
