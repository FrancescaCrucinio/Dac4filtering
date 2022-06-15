# distances
distances <- data.frame(aggregate(w1 ~ algo + runtime + N + run + d, data = df, FUN = "mean"))
distances$ks <- aggregate(ks ~ algo + runtime + N + run + d, data = df, FUN = "mean")$ks
time_means <- aggregate(runtime ~ algo + N + d, data = df, FUN= "mean" )
time_means <- time_means[order(time_means$algo, time_means$N, time_means$d), ]
distances <- distances[order(distances$algo, distances$N, distances$d), ]
distances$runtime_mean <- c(rep(time_means$runtime[1:11], each = 50), rep(time_means$runtime[12], each = 39), rep(time_means$runtime[13:26], each = 50))
distances_mean <- data.frame(aggregate(w1 ~ algo + N + d, data = distances, FUN = "mean"))
distances_mean$ks <- aggregate(ks ~ algo + N + d, data = distances, FUN = "mean")$ks
distances_mean <- merge(distances_mean, time_means, by=c("algo", "N", "d"))
cbPalette <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
# Wasserstein-1
ggplot(data = distances, aes(x = runtime_mean, y = w1, group = interaction(algo, N, d), fill = algo, colour = algo)) +
  geom_boxplot(coef = 6, width = 0.1, alpha = 0.1, lwd = 1) +
  # geom_point(data = distances_mean, shape = 4, lwd = 1, aes(x = runtime, y = w1, group = interaction(algo, N, d), fill = algo, colour = algo)) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~factor(d, levels = c("32", "256", "2048")), ncol = 1, nrow = 3) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  theme(axis.title.x=element_blank(),  axis.text = element_text(size=30),
        axis.title.y=element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30), legend.position="none", panel.spacing = unit(2, "lines"))
# ggsave("marginal_lgssm_w1.pdf", width = 20, height = 35, dpi = 300)
# Kolmogorov-Smirnov
ggplot(data = distances, aes(x = runtime_mean, y = ks, group = interaction(algo, N, d), fill = algo, colour = algo)) +
  geom_boxplot(coef = 10, width = 0.1, alpha = 0.1, lwd = 1) +
  # geom_point(data = distances_mean, shape = 4, lwd = 1, aes(x = runtime, y = ks, group = interaction(algo, N, d), fill = algo, colour = algo)) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~factor(d, levels = c("32", "256", "2048")), ncol = 1, nrow = 3) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  theme(axis.title.x=element_blank(), axis.text = element_text(size=30),
        axis.title.y=element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30), legend.position="none", panel.spacing = unit(2, "lines"))
# ggsave("marginal_lgssm_ks.pdf", width = 20, height = 35, dpi = 300)

p <- ggplot(data = distances, aes(x = runtime_mean, y = ks, group = interaction(algo, N, d), fill = algo, colour = algo)) +
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
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  theme(axis.title.x=element_blank(), axis.text = element_text(size=30),
        axis.title.y=element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30), legend.position="bottom", panel.spacing = unit(2, "lines"))
my_legend <- get_legend(p)
legend_p <- as_ggplot(my_legend)
# ggsave("lgssm_legend.pdf", width = 8, height = 1, dpi = 300)
