library(tidyr)
library(ggplot2)
set.seed(1234)
# get data
y <- read.csv(file="data/synthetic_data_msv_y")
true_x <- read.csv(file="data/synthetic_data_msv_x")

# dac results
Time.step <- 122
p <- 26
d <- 30
df <- read.csv("data/msv/means_dac_msv_syn_d30N500ID1", row.names = 1)
df$time <- 1:122
df$run <- 1
df$type <- rep(c("mean", "var", "min", "q1", "q2", "q3", "max"), each = Time.step)
# for (id in 2:50){
#   filename <- paste0("data/msv/means_dac_msv_syn_d30N500ID", id, sep = "")
#   dfnew <- read.csv(filename, row.names = 1)
#   dfnew$run <- id
#   dfnew$time <- 1:122
#   dfnew$type <- rep(c("mean", "var", "min", "q1", "q2", "q3", "max"), each = Time.step)
#   df <- rbind(df, dfnew)
# }
x_levels  <- vector(mode='character',length=30)
for (i in 1:30) {
  x_levels[i] <- paste0("X", i)
}
# df <- aggregate(. ~ time + type, data = df, FUN = "mean")
# df <- df[, c(3:32, 1:2, 33:34)]
dim_selection <- c(1, 20, 28)
df_plot <- gather(df[df$type == "mean", dim_selection], "dim")
colnames(df_plot)[2] <- "mean"
df_plot$var <-  gather(df[df$type == "var", dim_selection], "dim")[, 2]
df_plot$min <-  gather(df[df$type == "min", dim_selection], "dim")[, 2]
df_plot$q1 <-  gather(df[df$type == "q1", dim_selection], "dim")[, 2]
df_plot$q2 <-  gather(df[df$type == "q2", dim_selection], "dim")[, 2]
df_plot$q3 <-  gather(df[df$type == "q3", dim_selection], "dim")[, 2]
df_plot$max <-  gather(df[df$type == "max", dim_selection], "dim")[, 2]
df_plot$time <- 1:122
df_plot_true_x <- gather(true_x[2:123, dim_selection], "dim")
df_plot_true_x$time <- 1:122
df_plot_true_x$N <- "truth"
colnames(y) <- x_levels[1:p]
df_plot_y <- gather(y[, dim_selection[dim_selection %in% 1:26]], "dim")
df_plot_y$time <- 1:122
df_plot_y$N <- "obs"

ggplot(data = df_plot, aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), fill = "grey70") +
  geom_line(color = "black") +
  geom_line(data = df_plot_true_x, aes(x = time, y = value), color = "blue") +
  geom_line(data = df_plot_y, aes(x = time, y = value), color = "red") +
  facet_wrap(~factor(dim, levels = x_levels[dim_selection]), scales = "free", nrow = 3) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_text(size=5),
        strip.text.x = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=5))


# for (id in 1:50){
#   filename <- paste0("data/msv/means_dac_msv_syn_d30N500ID", id, sep = "")
#   dfnew <- read.csv(filename, row.names = 1)
#   dfnew$N <- "5*10^2"
#   dfnew$run <- id
#   dfnew$algo <- "dac"
#   dfnew$time <- 1:122
#   dfnew$type <- rep(c("mean", "var", "min", "q1", "q2", "q3", "max"), each = Time.step)
#   df <- rbind(df, dfnew)
# }
# for (id in 1:50){
#   filename <- paste0("data/msv/means_dac_msv_syn_d30N1000ID", id, sep = "")
#   dfnew <- read.csv(filename, row.names = 1)
#   dfnew$N <- "10^3"
#   dfnew$run <- id
#   dfnew$algo <- "dac"
#   dfnew$time <- 1:122
#   df <- rbind(df, dfnew)
# }
#
#
# tmp <- aggregate(. ~ algo + N + time, data = df, FUN = "mean")
# tmp <- tmp[order(tmp$algo, tmp$N), ]
# x_levels  <- vector(mode='character',length=30)
# for (i in 1:30) {
#   x_levels[i] <- paste0("X", i)
# }
# df_plot <- gather(tmp[, 4:33], "dim")
# df_plot$time <- 1:122
# df_plot$N <- rep(c("10^2", "5*10^2", "10^3"), each = 122)
# df_plot_true_x <- gather(true_x[2:123, ], "dim")
# df_plot_true_x$time <- 1:122
# df_plot_true_x$N <- "truth"
# colnames(y) <- x_levels[1:p]
# df_plot_y <- gather(y, "dim")
# df_plot_y$time <- 1:122
# df_plot_y$N <- "obs"
# ggplot(data = df_plot, aes(x = time, y = value, group = N, color = N)) +
#   geom_line(aes(group = N)) +
#   geom_line(data = df_plot_true_x, aes(x = time, y = value), color = "black") +
#   geom_line(data = df_plot_y, aes(x = time, y = value), color = "gray") +
#   facet_wrap(~factor(dim, levels = x_levels), nrow = 15, ncol = 2, scales = "free", dir = "v") +
#   theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
#         axis.text = element_text(size=5),
#         strip.text.x = element_blank(),
#         legend.title = element_blank(), legend.text=element_text(size=20),
#         text = element_text(size=5))


tmp[tmp$type == "mean", 1]
tmp[tmp$type == "min", 1]
