# Get stats to compare resampling schemes for linear Gaussian SSM
library(ggpubr)
d <- 128
# read data
df <- rbind(read.csv(paste0("data/resampling_tempering/light_marginal_resampling_comparison_d", d, "N100ID1timeinterval1")),
            read.csv(paste0("data/resampling_tempering/ada_light_marginal_resampling_comparison_d", d, "N100ID1timeinterval1")),
            read.csv(paste0("data/resampling_tempering/mix_marginal_resampling_comparison_d", d, "N100ID1timeinterval1")),
            read.csv(paste0("data/resampling_tempering/lc_marginal_resampling_comparison_d", d, "N100ID1timeinterval1")))
df$N <- "10^2"
for (id in 2:50){
  dfnew <- rbind(read.csv(paste0("data/resampling_tempering/light_marginal_resampling_comparison_d", d, "N100ID", id, "timeinterval1", sep = "")),
                 read.csv(paste0("data/resampling_tempering/ada_light_marginal_resampling_comparison_d", d, "N100ID", id, "timeinterval1", sep = "")),
                 read.csv(paste0("data/resampling_tempering/mix_marginal_resampling_comparison_d", d, "N100ID", id, "timeinterval1", sep = "")),
                 read.csv(paste0("data/resampling_tempering/lc_marginal_resampling_comparison_d", d, "N100ID", id, "timeinterval1", sep = "")))
  dfnew$N <- "10^2"
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  dfnew <- rbind(read.csv(paste0("data/resampling_tempering/light_marginal_resampling_comparison_d", d, "N1000ID", id, "timeinterval1", sep = "")),
                 read.csv(paste0("data/resampling_tempering/ada_light_marginal_resampling_comparison_d", d, "N1000ID", id, "timeinterval1", sep = "")),
                 #read.csv(paste0("data/resampling_tempering/mix_marginal_resampling_comparison_d", d, "N1000ID", id, "timeinterval1", sep = "")),
                 read.csv(paste0("data/resampling_tempering/lc_marginal_resampling_comparison_d", d, "N1000ID", id, "timeinterval1", sep = "")))
  dfnew$N <- "10^3"
  df <- rbind(df, dfnew)
}
for (id in 1:50){
  dfnew <- rbind(read.csv(paste0("data/resampling_tempering/light_marginal_resampling_comparison_d", d, "N10000ID", id, "timeinterval1", sep = "")),
                 read.csv(paste0("data/resampling_tempering/ada_light_marginal_resampling_comparison_d", d, "N10000ID", id, "timeinterval1", sep = "")),
                 read.csv(paste0("data/resampling_tempering/mix_marginal_resampling_comparison_d", d, "N10000ID", id, "timeinterval1", sep = "")),
                 read.csv(paste0("data/resampling_tempering/lc_marginal_resampling_comparison_d", d, "N10000ID", id, "timeinterval1", sep = "")))

  dfnew$N <- "10^4"
  df <- rbind(df, dfnew)
}
df <- df[, -1]
d <- ncol(df) - 3
df$d_means <- rowMeans(df[, 1:d])
df <- df[order(df$algo, df$N),]
time_means <- aggregate(elapsed ~ algo + N, data = df, FUN = "mean" )
time_means <- time_means[order(time_means$algo, time_means$N), ]
df$runtime <- rep(time_means$elapsed, each = 50)

data_df <- df[, (d+2):ncol(df)]
data_df$d <- d
write.csv(data_df, file = paste0("data/resampling_tempering/marginal_resampling_comparison_d", d, sep = ""))
