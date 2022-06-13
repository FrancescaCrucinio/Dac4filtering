### Spatial model stats
library(ggplot2)
# dimension
d <- 8
# parameters
sigmaX <- 1
nu <- 10
tau_diag <- 1
tau <- -0.25
df <- data.frame()

for (Nparticles in c(100, 500, 1000)) {
  df_dac <- read.csv(paste0("data/spatial/same_nt_corrected_stats_dac_spatial_tau", -tau, "d", d, "N", Nparticles),
                     row.names = 1)
  df <- rbind(df, df_dac)
}

dim <- 1
t <- 10
df_plot <- df[df$t == t & df$dim == dim,]
ggplot(data = df_plot, aes(x=runtime, y=var, group = N))+
  geom_boxplot(coef = 10, width = 0.1, alpha = 0.1, lwd = 1) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30))


ground_truth <- unname(data.matrix(read.csv(paste0("data/spatial/corrected_data_truth_spatial_tau_", -tau, "d", d, "ID", 1),
                                            row.names = 1, nrows=d, skip=(Time.step)*d)))
df_mse <- df[df$t == 10, ]
df_mse$col <- ceiling(df_mse$dim/d)
df_mse$row <- df_mse$dim - (df_mse$col - 1)*d
df_mse$mse <- (df_mse$mean - rep(c(ground_truth), times = 3*50))^2
tmp <- aggregate(mse ~  N + col + row, data = df_mse, FUN = "mean")
ggplot(data = tmp, aes(x = col, y = row, fill = mse)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~N, ncol = 3, nrow = 1) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.text.x = element_blank(), panel.spacing = unit(1, "lines"),
        legend.title = element_blank(), legend.text=element_text(size=25),
        text = element_text(size=20), aspect.ratio = 1)
