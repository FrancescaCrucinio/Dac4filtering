### Spatial model stats
library(ggplot2)
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
  df_dac$type <- "dac"
  df <- rbind(df, df_dac)
}
df_bpf <- read.csv(paste0("data/spatial/new_stats_bpf_spatial_tau", -tau, "d", d, "N", 100000),
                   row.names = 1)
dim <- 3
t <- 10
df_plot <- df[df$t == t & df$dim == dim,]
df_bpf_plot <- df_bpf[df_bpf$t == t & df_bpf$dim == dim,]
ggplot(data = df_plot, aes(x=N, y=mean, group = N))+
  geom_boxplot(coef = 10, width = 0.1, alpha = 0.3, lwd = 1, color = "dark gray", fill = "dark gray") +
  geom_hline(yintercept = mean(df_bpf_plot$mean), lwd = 1, linetype = "dashed") +
  geom_hline(yintercept = quantile(df_bpf_plot$mean, 0.25), lwd = 1, linetype = "dashed") +
  geom_hline(yintercept = quantile(df_bpf_plot$mean, 0.75), lwd = 1, linetype = "dashed") +
    scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30))
# ggsave("spatial2_boxplot_node12_mean.pdf", width = 8.5, height = 6, dpi = 300)
