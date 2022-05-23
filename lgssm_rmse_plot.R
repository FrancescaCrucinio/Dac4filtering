# RMSE for linear Gaussian SSM
library(ggpubr)
# read data
df <- read.csv(paste0("data/lgssm_tempering/lgssm_d32N100ID1"))
df$N <- "10^2"
df$d <- "32"
df$run <- 1
dfnew <- read.csv(paste0("data/lgssm_tempering/marginal_lgssm_d32N100ID1"))
dfnew$N <- "10^2"
dfnew$d <- "32"
dfnew$run <- 1
df <- rbind(df, dfnew)
dfnew <- read.csv(paste0("data/lgssm_tempering/lgssm_d256N100ID1"))
dfnew$N <- "10^2"
dfnew$d <- "256"
dfnew$run <- 1
df <- rbind(df, dfnew)
dfnew <- read.csv(paste0("data/lgssm_tempering/marginal_lgssm_d256N100ID1"))
dfnew$N <- "10^2"
dfnew$d <- "256"
dfnew$run <- 1
df <- rbind(df, dfnew)
dfnew <- read.csv(paste0("data/lgssm_tempering/lgssm_d2048N100ID1"))
dfnew$N <- "10^2"
dfnew$d <- "2048"
dfnew$run <- 1
df <- rbind(df, dfnew)
dfnew <- read.csv(paste0("data/lgssm_tempering/marginal_lgssm_d2048N100ID1"))
dfnew$N <- "10^2"
dfnew$d <- "2048"
dfnew$run <- 1
df <- rbind(df, dfnew)
for (id in 2:50){
  filename <- paste0("data/lgssm_tempering/lgssm_d32N100ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^2"
  dfnew$d <- "32"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  dfnew <- read.csv(paste0("data/lgssm_tempering/marginal_lgssm_d32N100ID", id, sep = ""))
  dfnew$N <- "10^2"
  dfnew$d <- "32"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/lgssm_d256N100ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^2"
  dfnew$d <- "256"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  dfnew <- read.csv(paste0("data/lgssm_tempering/marginal_lgssm_d256N100ID", id, sep = ""))
  dfnew$N <- "10^2"
  dfnew$d <- "256"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/lgssm_d2048N100ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^2"
  dfnew$d <- "2048"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/marginal_lgssm_d2048N100ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^2"
  dfnew$d <- "2048"
  dfnew$run <- id
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  filename <- paste0("data/lgssm_tempering/lgssm_d32N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^3"
  dfnew$d <- "32"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/marginal_lgssm_d32N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^3"
  dfnew$d <- "32"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/lgssm_d256N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^3"
  dfnew$d <- "256"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/marginal_lgssm_d256N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^3"
  dfnew$d <- "256"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/lgssm_d2048N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^3"
  dfnew$d <- "2048"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/marginal_lgssm_d2048N1000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^3"
  dfnew$d <- "2048"
  dfnew$run <- id
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  filename <- paste0("data/lgssm_tempering/lgssm_d32N10000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^4"
  dfnew$d <- "32"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  dfnew <- read.csv(paste0("data/lgssm_tempering/marginal_lgssm_d32N10000ID", id, sep = ""))
  dfnew$N <- "10^4"
  dfnew$d <- "32"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/lgssm_d256N10000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^4"
  dfnew$d <- "256"
  dfnew$run <- id
  df <- rbind(df, dfnew)
  filename <- paste0("data/lgssm_tempering/marginal_lgssm_d256N10000ID", id, sep = "")
  dfnew <- read.csv(filename)
  dfnew$N <- "10^4"
  dfnew$d <- "256"
  dfnew$run <- id
  df <- rbind(df, dfnew)
}
df <- df[, -1]
df <- df[df$algo != "dac_ada", ]
Time.step <- ncol(df) - 7
colnames(df)[(Time.step+1):(Time.step+3)] <- c("w1", "ks", "runtime")

tmp <- aggregate(. ~ algo + N + d, data = df, FUN = "mean")
rmse_data <- data.frame(rep(1:Time.step, times = 29), rep(tmp$algo, each = 100), rep(tmp$N, each = 100), rep(tmp$d, each = 100))
colnames(rmse_data) <- c("Time.step", "algo", "N", "d")
rmse_data <- rmse_data[order(rmse_data$algo, rmse_data$N, rmse_data$d), ]
tmp <- tmp[order(tmp$algo, tmp$N, tmp$d), ]
rmse_data$rmse <- as.vector(t(as.matrix(tmp[, 4:103])))
ggplot(data = rmse_data, aes(x = Time.step, y = rmse, group = algo, colour = algo)) +
  geom_line(size = 2) +
  scale_y_continuous(trans='log10') +
  facet_wrap(~factor(interaction(N, d), levels=c("10^2.32", "10^3.32", "10^4.32", "10^2.256", "10^3.256", "10^4.256", "10^2.2048", "10^3.2048")), ncol = 3, nrow = 3) +
  scale_color_manual(values=c("#F8766D", "#999999", "#00BA38", "#619CFF")) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_text(size=20), strip.text.x = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15))
# ggsave("marginal_lgssm_rmse.pdf", width = 12, height = 8, dpi = 300)
