output_stats <- function(ID, d, Time.step, Nparticles, M, model){
  # get exact means, variances and marginals
  if(model == "lgssm"){
    filename <- paste0("data/data_lgssm_d", d, "ID", ID)
  }
  data <- read.csv(filename, row.names = 1)
  true_means <- data.matrix(data[data$type == "mean", 1:d])
  true_variances <- data.matrix(data[data$type == "var", 1:d])
  marginals <- matrix(0, nrow = 10^5, ncol = d)
  for(i in 1:d){
    marginals[, i] <- rnorm(10^5, mean = true_means[Time.step, i], sd = sqrt(true_variances[Time.step, i]))
  }

  # number of timeintervals (aka data sets)
  ntimeinterval <- Time.step/10
  data_dac <- data.frame()
  data_nsmc <- data.frame()
  data_stpf <- data.frame()
  for (i in 1:ntimeinterval) {
    if(model == "lgssm"){
      filename <- paste0("data/dac_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", i)
      data_dac <- rbind(data_dac, read.csv(filename, row.names = 1))
      filename <- paste0("data/nsmc_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", i)
      data_nsmc <- rbind(data_nsmc, read.csv(filename, row.names = 1))
      filename <- paste0("data/stpf_lgssm_d", d, "N", Nparticles, "ID", ID, "timeinterval", i)
      data_stpf <- rbind(data_stpf, read.csv(filename, row.names = 1))
    }
  }
  runtime_dac <- sum(data_dac$runtime[seq(from = 1, to = Time.step*Nparticles, by = Nparticles)])
  runtime_nsmc <- sum(data_nsmc$runtime[seq(from = 1, to = Time.step*Nparticles, by = Nparticles)])
  runtime_stpf <- sum(data_stpf$runtime[seq(from = 1, to = Time.step*Nparticles*M, by = M*Nparticles)])
  res_dac <- unname(data.matrix(data_dac[data_dac$timestep == Time.step, 1:d]))
  res_nsmc <- unname(data.matrix(data_nsmc[data_nsmc$timestep == Time.step, 1:d]))
  res_stpf <- array(c(data.matrix(data_stpf[data_stpf$timestep == Time.step, 1:d])), dim = c(Nparticles, M, d))
  # distances (last time step only)
  ks_dac <- apply(rbind(res_dac, marginals), ks_dist, N = Nparticles, MARGIN = 2)
  w1_dac <- apply(rbind(res_dac, marginals), w1_dist, N = Nparticles, MARGIN = 2)
  ks_nsmc <- apply(rbind(res_nsmc, marginals), ks_dist, N = Nparticles, MARGIN = 2)
  w1_nsmc <- apply(rbind(res_nsmc, marginals), w1_dist, N = Nparticles, MARGIN = 2)
  ks_stpf <- apply(rbind(matrix(res_stpf, ncol = d, nrow = Nparticles*M), marginals), ks_dist, N = M*Nparticles, MARGIN = 2)
  w1_stpf <- apply(rbind(matrix(res_stpf, ncol = d, nrow = Nparticles*M), marginals), w1_dist, N = M*Nparticles, MARGIN = 2)
  df_distances <- data.frame(cbind(c(ks_dac, ks_nsmc, ks_stpf), c(w1_dac, w1_nsmc, w1_stpf)))
  df_distances$algo <- rep(c("dac", "nsmc", "stpf"), each = d)
  df_distances$runtime <- rep(c(runtime_dac, runtime_nsmc, runtime_stpf), each = d)
  colnames(df_distances)[1:2] <- c("ks", "w1")
  if(model == "lgssm"){
    filename <- paste0("/storage/u1693998/out/distances_lgssm_d", d, "N", Nparticles, "ID", ID)
  }
  write.csv(x=df_distances, file=filename)
  # mean and variance
  m_dac <- matrix(0, nrow = d, ncol = Time.step)
  v_dac <- matrix(0, nrow = d, ncol = Time.step)
  m_nsmc <- matrix(0, nrow = d, ncol = Time.step)
  v_nsmc <- matrix(0, nrow = d, ncol = Time.step)
  m_stpf <- matrix(0, nrow = d, ncol = Time.step)
  v_stpf <- matrix(0, nrow = d, ncol = Time.step)
  for (t in 1:Time.step) {
    res_dac <- data.matrix(data_dac[data_dac$timestep == t, 1:d])
    res_nsmc <- data.matrix(data_nsmc[data_nsmc$timestep == t, 1:d])
    res_stpf <- array(c(data.matrix(data_stpf[data_stpf$timestep == t, 1:d])), dim = c(Nparticles, M, d))
    m_dac[, t] <- (colMeans(res_dac) - true_means[t, ])^2/true_variances[t, ]
    v_dac[, t] <- (colVars(res_dac) - true_variances[t, ])^2
    m_nsmc[, t] <- (colMeans(res_nsmc) - true_means[t, ])^2/true_variances[t, ]
    v_nsmc[, t] <- (colVars(res_nsmc) - true_variances[t, ])^2
    m_stpf[, t] <- (colMeans(res_stpf, dims = 2) - true_means[t, ])^2/true_variances[t, ]
    v_stpf[, t] <- (colMeans(res_stpf^2, dims = 2) - m_stpf[, t]^2 - true_variances[t, ])^2
  }
  if(model == "lgssm"){
    filename <- paste0("/storage/u1693998/out/rmse_lgssm_d", d, "N", Nparticles, "ID", ID)
  }
  df_rmse <- data.frame(rbind(m_dac, m_nsmc, m_stpf, v_dac, v_nsmc, v_stpf))
  df_rmse$algo <- rep(c("dac", "nsmc", "stpf"), each = 2*d)
  write.csv(x=df_rmse, file=filename)
}
# component <- 25
# library(ks)
# kde_dac <- kde(x = res_dac[, component])
# kde_nsmc <- kde(x = res_nsmc[, component], eval.points = kde_dac$eval.points)
# kde_stpf <- kde(x = c(res_stpf[, , component]), eval.points = kde_dac$eval.points)
# plot(kde_dac$eval.points, dnorm(kde_dac$eval.points, mean = true_means[Time.step, component],
#                                 sd = sqrt(true_variances[Time.step, component])), type = "l", ylim = c(0, 1.2))
# lines(kde_dac$eval.points, kde_dac$estimate, type = "l", col = "blue")
# lines(kde_dac$eval.points, kde_nsmc$estimate, type = "l", col = "red")
# lines(kde_dac$eval.points, kde_stpf$estimate, type = "l", col = "green")
#
# plot(1:32, df_distances$ks[1:32], type = "l", col = "blue", ylim = c(0, 0.5))
# lines(1:32, df_distances$ks[33:64], type = "l", col = "red")
# lines(1:32, df_distances$ks[65:96], type = "l", col = "green")
