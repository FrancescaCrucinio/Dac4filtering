library(tidyr)
# get parameters
load("~/Documents/Dac4filtering/data/msv_parameters.RData")
# get true data
true_x <- read.csv(file="data/synthetic_data_msv_x")
# read dac data
df <- read.csv("data/msv/means_dac_msv_syn_d30N100ID1", row.names = 1)
df$N <- "10^2"
df$run <- 1
df$time <- 1:122
df$algo <- "dac"
for (id in 2:50){
  filename <- paste0("data/msv/means_dac_msv_syn_d30N100ID", id, sep = "")
  dfnew <- read.csv(filename, row.names = 1)
  dfnew$N <- "10^2"
  dfnew$run <- id
  dfnew$algo <- "dac"
  dfnew$time <- 1:122
  df <- rbind(df, dfnew)
}
tmp <- aggregate(. ~ algo + N + time, data = df, FUN = "mean")
# df_plot <- gather(rbind(true_x[2:123,], tmp[, 4:33]), "dim")
# df_plot$time <- 1:122
# df_plot$type <- rep(c("truth", "filter"), each = 122)
# x_levels  <- vector(mode='character',length=30)
# for (i in 1:30) {
#   x_levels[i] <- paste0("X", i)
# }
# ggplot(data = df_plot, aes(x = time, y = value, group = type, color = type)) +
#   geom_line(aes(group = type)) +
#   facet_wrap(~factor(dim, levels = x_levels), nrow = 15, ncol = 2, scales = "free", dir = "v") +
#   theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
#         axis.text = element_text(size=5),
#         # strip.text.x = element_blank(),
#         legend.title = element_blank(), legend.text=element_text(size=20),
#         text = element_text(size=5))


# run factorstochvol
res <- fsvsample(y, factors = 4, draws = 50000, zeromean = TRUE, thin = 10, quiet = TRUE, restric = "auto")
