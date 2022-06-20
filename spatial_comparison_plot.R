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
  df_dac <- read.csv(paste0("data/spatial/same_nt_corrected_stats_dac_spatial_tau", -tau, "d", d, "N", Nparticles),
                     row.names = 1)
  df_dac$type <- "dac"
  df_bpf <- read.csv(paste0("data/spatial/same_nt_corrected_stats_bpf_spatial_tau", -tau, "d", d, "N", Nparticles),
                     row.names = 1)
  df_bpf$type <- "bpf"
  df <- rbind(df, df_dac, df_bpf)
}

dim <- 2
t <- 10
df_plot <- df[df$t == t & df$dim == dim,]
ggplot(data = df_plot, aes(x=N, y=mean, group = interaction(N, type), color = type, fill = type))+
  geom_boxplot(coef = 10, width = 0.1, alpha = 0.1, lwd = 1) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30))
# ggsave("spatial2_boxplot_node22_mean.pdf", width = 8.5, height = 6, dpi = 300)
