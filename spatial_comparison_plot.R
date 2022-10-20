### Spatial model stats
library(ggplot2)
library(ggpubr)
# dimension
d <- 2
# parameters
sigmaX <- 1
nu <- 10
tau_diag <- 1
tau <- -0.25
Time.step <- 10
df <- data.frame()

for (Nparticles in c(100, 500, 1000, 5000, 10000)) {
  df_dac <- read.csv(paste0("data/spatial/new_stats_dac_spatial_tau", -tau, "d", d, "N", Nparticles),
                     row.names = 1)
  df_dac$type <- "dac-ada"
  df <- rbind(df, df_dac)
}
for (Nparticles in c(100, 500, 1000, 5000)) {
  df_dac <- read.csv(paste0("data/spatial/nonadaptive_stats_dac_spatial_tau", -tau, "d", d, "N", Nparticles),
                     row.names = 1)
  df_dac$type <- "dac"
  df <- rbind(df, df_dac)
}
df_bpf <- read.csv(paste0("data/spatial/new_stats_bpf_spatial_tau", -tau, "d", d, "N", 100000),
                   row.names = 1)
dim <- 4
t <- 10
df_plot <- df[df$t == t & df$dim == dim,]
df_bpf_plot <- df_bpf[df_bpf$t == t & df_bpf$dim == dim,]
ggplot(data = df_plot, aes(x=N, y=mean, group = interaction(N, type), color = type, fill = type))+
  geom_boxplot(coef = 10, width = 0.1, alpha = 0.3, lwd = 1) +
  geom_hline(yintercept = mean(df_bpf_plot$mean), lwd = 1, linetype = "dashed") +
  geom_hline(yintercept = quantile(df_bpf_plot$mean, 0.25), lwd = 1, linetype = "dashed") +
  geom_hline(yintercept = quantile(df_bpf_plot$mean, 0.75), lwd = 1, linetype = "dashed") +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30), legend.position="none")
# ggsave("bpf_spatial2_boxplot_node22_mean.pdf", width = 8.5, height = 6, dpi = 300)

p <- ggplot(data = df_plot, aes(x=N, y=mean, group = interaction(N, type), color = type, fill = type))+
  geom_boxplot(coef = 10, width = 0.1, alpha = 0.3, lwd = 1) +
  geom_hline(yintercept = mean(df_bpf_plot$mean), lwd = 1, linetype = "dashed") +
  geom_hline(yintercept = quantile(df_bpf_plot$mean, 0.25), lwd = 1, linetype = "dashed") +
  geom_hline(yintercept = quantile(df_bpf_plot$mean, 0.75), lwd = 1, linetype = "dashed") +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30), legend.position="bottom")
my_legend <- get_legend(p)
legend_p <- as_ggplot(my_legend)
# ggsave("spatial_legend_bpf.pdf", width = 8, height = 1, dpi = 300)


df_final <- df[df$t == t, c(1, 7:11)]
runtime_means <- aggregate(runtime ~ type + N, data = df_final, FUN = mean)
ggplot(data = runtime_means, aes(x=N, y= runtime, group = type, color = type)) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_point(lwd = 4, aes(shape = as.factor(N))) +
  geom_line(lwd = 1) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_text(size=20), strip.text.x = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15))
# ggsave("spatial2_cost_per_N.pdf", width = 8.5, height = 6, dpi = 300)
