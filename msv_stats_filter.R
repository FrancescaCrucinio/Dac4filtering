library(tidyr)
library(ggplot2)
# get data
y <- read.csv(file="data/synthetic_data_msv_y")
true_x <- read.csv(file="data/synthetic_data_msv_x")
Time.step <- dim(true_x)[1]
d <- dim(true_x)[2]
# get stpf results
filename <- "data/msv/means_stpf_msv_syn_d128N50ID1"
df <- read.csv(filename, row.names = 1)
# get dac results
filename <- "data/msv/means_dac_msv_syn_d128N100ID1"
df_dac <- read.csv(filename, row.names = 1)
# plot
dim_selection <- c(1, 20, 28)
df_plot <- rbind(gather(df[df$type == "mean", dim_selection], "dim"), gather(df_dac[df_dac$type == "mean", dim_selection], "dim"))
df_plot$time <- 1:Time.step
df_plot$q1 <-  c(gather(df[df$type == "q1", dim_selection], "dim")[, 2], gather(df_dac[df_dac$type == "q1", dim_selection], "dim")[, 2])
df_plot$q3 <-  c(gather(df[df$type == "q3", dim_selection], "dim")[, 2], gather(df_dac[df_dac$type == "q3", dim_selection], "dim")[, 2])
df_plot$algo <- rep(c("stpf", "dac"), each = Time.step*length(dim_selection))
df_plot$truth <- gather(true_x[, dim_selection], "dim")[, 2]
x_levels  <- vector(mode='character',length = d)
for (i in 1:d) {
  x_levels[i] <- paste0("X", i)
}
cbPalette <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
ggplot(data = df_plot, aes(x = time, y = value, group = algo, color = algo)) +
  # geom_line(aes(x = time, y = q1, group = algo, color = algo), linetype = "dashed") +
  # geom_line(aes(x = time, y = q3, group = algo, color = algo), linetype = "dashed") +
  geom_line(aes(color = algo), linetype = 2) +
  geom_line(aes(x = time, y = truth, color = "truth"), lwd = 1) +
  facet_wrap(~factor(dim, levels = x_levels[dim_selection]), scales = "free", nrow = 3) +
  scale_color_manual(values=c(cbPalette[1], cbPalette[5], "black")) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_text(size=15),
        strip.text.x = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15))
# ggsave("marginal_msv_synthetic_filtering.pdf", width = 12, height = 8, dpi = 300)

dim_selection <- 1
df_plot2 <- data.frame(c(df[df$type == "mean", dim_selection], df_dac[df_dac$type == "mean", dim_selection]))
colnames(df_plot2) <- "mean"
df_plot2$time <- 1:Time.step
df_plot2$q1 <-  c(df[df$type == "q1", dim_selection], df_dac[df_dac$type == "q1", dim_selection])
df_plot2$q3 <-  c(df[df$type == "q3", dim_selection], df_dac[df_dac$type == "q3", dim_selection])
df_plot2$iqr <-  c(df[df$type == "q3", dim_selection] - df[df$type == "q1", dim_selection],
                   df_dac[df_dac$type == "q3", dim_selection] - df_dac[df_dac$type == "q1", dim_selection])
df_plot2$algo <- rep(c("stpf", "dac"), each = Time.step*length(dim_selection))
df_plot2$truth <- true_x[, dim_selection]
ggplot(data = df_plot2, aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = mean - 1.5*iqr, ymax = mean + 1.5*iqr), fill = "grey70") +
  geom_line() +
  geom_line(aes(x = time, y = truth, color = "truth"), lwd = 1) +
  facet_wrap(~factor(algo), nrow = 2) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_text(size=15), strip.text = element_text(size = 20),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15))
# ggsave("marginal_msv_synthetic_ribbon.pdf", width = 12, height = 8, dpi = 300)
