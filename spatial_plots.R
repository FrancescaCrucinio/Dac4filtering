### Spatial model stats
library(ggplot2)
# dimension
d <- 16
# parameters
sigmaX <- 1
nu <- 10
tau_diag <- 1
tau <- -0.25
Time.step <- 10
df <- data.frame()

for (Nparticles in c(100, 500, 1000)) {
  df_dac <- read.csv(paste0("data/spatial/adaptive_stats_dac_spatial_tau", -tau, "d", d, "N", Nparticles),
                     row.names = 1)
  df <- rbind(df, df_dac)
}
# last time step
df_end <- df[df$t == Time.step, ]
df_iqr_mean <- aggregate(mean ~ dim + N, data = df_end, FUN = IQR)

ggplot(data = df_iqr_mean, aes(x=N, y=mean, group=dim, color=dim)) +
  geom_line(alpha=0.5) +
  geom_abline(slope = -0.5, intercept = 10^(0.01), size = 2, linetype = "dashed") +
  # geom_smooth(method='lm', formula= y~x, fill = NA, aes(group=1)) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  labs(color='d') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_text(size=20), strip.text.x = element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=30))
# ggsave("spatial8_iqr.pdf", width = 8.5, height = 4, dpi = 300)


dim <- 10
t <- 10
df_plot <- df[df$t == t & df$dim == dim,]
ggplot(data = df_plot, aes(x=runtime, y=mean, group = N))+
  geom_boxplot(coef = 10, alpha = 0.1, lwd = 1) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30))
# ggsave("spatial8_boxplot_node22_mean.pdf", width = 8.5, height = 6, dpi = 300)
