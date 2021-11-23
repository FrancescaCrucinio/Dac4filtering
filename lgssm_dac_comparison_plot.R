
mse <- rep(0, times = 3*Nrep)
rmse <- rep(0, times = 3*Nrep)
for (j in 1:Nrep) {
  mse[j] <- rowMeans(df_dac[j*Time.step, 1:d])
  rmse[j] <- rowMeans(df_dac[j*Time.step, 1:d]/true_variances)
  mse[Nrep+j] <- rowMeans(df_dac_mix[j*Time.step, 1:d])
  rmse[Nrep+j] <- rowMeans(df_dac_mix[j*Time.step, 1:d]/true_variances)
  mse[2*Nrep+j] <- rowMeans(df_dac_light[j*Time.step, 1:d])
  rmse[2*Nrep+j] <- rowMeans(df_dac_light[j*Time.step, 1:d]/true_variances)
}

g <- as.factor(rep(c("dac", "mix", "light"), each = Nrep))
times <- rep(res$total_time/res$n_itr, each = Nrep)
memory <- rep(res$mem_alloc/1e6, each = Nrep)
# time
df <- data.frame(x1 = times, x2 = memory, y = mse, g)
ggplot(data = df, aes(x = x1, y = y, group = g, fill = g)) +
  geom_boxplot(aes(x = x1, y= y), coef = 6) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=15))
# memory
ggplot(data = df, aes(x = x2, y = y, group = g, fill = g)) +
  geom_boxplot(aes(x = x2, y= y), coef = 6) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=15))
