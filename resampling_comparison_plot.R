# Plot of MSE for different resampling schemes for linear Gaussian SSM
library(ggpubr)
# read data
df <- rbind(read.csv(paste0("data/resampling_tempering/resampling_comparison_d8")),
            read.csv(paste0("data/resampling_tempering/resampling_comparison_d128")),
            read.csv(paste0("data/resampling_tempering/resampling_comparison_d1024")))
# time
# p <-
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
  facet_wrap(~factor(d, levels = c("8", "128", "1024")), ncol = 1, nrow = 3) +
  theme(axis.title.x=element_blank(), axis.text = element_text(size=30),
        axis.title.y=element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=30),
        text = element_text(size=30), legend.position="none", panel.spacing = unit(2, "lines"))
# my_legend <- get_legend(p)
# as_ggplot(my_legend)
# ggsave("resampling_legend.pdf", width = 6, height = 1, dpi = 300)
# ggsave("resampling_comparison.pdf", width = 40, height = 40, dpi = 300)

