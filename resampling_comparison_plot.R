# Plot of MSE for different resampling schemes for linear Gaussian SSM
library(ggpubr)
# read data
df <- rbind(read.csv(paste0("data/resampling_tempering/marginal_resampling_comparison_d128")))
ggplot(data = df, aes(x = runtime, y = d_means, group = interaction(algo, N), fill = algo, colour = algo)) +
  geom_boxplot(coef = 6, width = 5) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.title.x=element_blank(), axis.text = element_text(size=40),
        axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=100),
        legend.position = "top", legend.key.width = unit(5, 'cm'))
# ggsave("marginal_resampling_comparison.pdf", width = 40, height = 20, dpi = 300)
