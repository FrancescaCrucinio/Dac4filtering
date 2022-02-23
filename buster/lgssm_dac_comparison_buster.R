devtools::load_all("/storage/u1693998/Dac4filtering")
### Linear Gaussian SSM -- comparison of dac and dac with mixture reweighting (both lightweight and full cost)
ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
set.seed(1234*ID)
# dimension
d <- 1024
# initial state
mu0 <- rep(0, times = d)
Sigma0 <- diag(x = 1, d, d)
# parameters
tau <- 1
lambda <- 1
sigmaY <- 0.5^2
Time.step <- 10

filename <- paste0("/storage/u1693998/data/data_lgssm_d", d, "ID", ID)
# filename <- paste0("data/data_lgssm_d", d, "ID", ID)
data <- read.csv(filename, row.names = 1)
y <- data.matrix(data[1:Time.step, 1:d], rownames.force = NA)
true_means <- data.matrix(data[data$type == "mean", 1:d])
true_variances <- data.matrix(data[data$type == "var", 1:d])
### DAC
Nparticles <- 100
df <- data.frame()

# NO CROSSOVER
# dac
x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
tic()
res_dac <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "lc")
runtime <- toc()
se <- (res_dac$m[Time.step, ] - true_means[Time.step, ])^2
df <- data.frame(rbind(df, t(c(se, runtime))))
# mix
# tic()
# res_dac_mix <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "mix")
# runtime <- toc()
# se <- (res_dac_mix$m[Time.step, ] - true_means[Time.step, ])^2
# df <- data.frame(rbind(df, t(c(se, runtime))))
# light
tic()
res_dac_light <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "light")
runtime <- toc()
se <- (res_dac_light$m[Time.step, ] - true_means[Time.step, ])^2
df <- data.frame(rbind(df, t(c(se, runtime))))
# adaptive light
tic()
res_dac_light_v2 <- dac_time_lgssm(tau, lambda, sigmaY, Nparticles, x0, y, method = "ada")
runtime <- toc()
se <- (res_dac_light_v2$m[Time.step, ] - true_means[Time.step, ])^2
mean(se)
df <- data.frame(rbind(df, t(c(se, runtime))))
# # CROSSOVER
# # dac
# x0 <- mvrnorm(n = Nparticles, mu0, Sigma0)
# tic()
# res_dac <- dac_time_lgssm_crossover(tau, lambda, sigmaY, Nparticles, x0, y, method = "lc")
# runtime <- toc()
# se <- (res_dac$m[Time.step, ] - true_means[Time.step, ])^2
# df <- data.frame(rbind(df, t(c(se, runtime))))
# # mix
# # tic()
# # res_dac_mix <- dac_time_lgssm_crossover(tau, lambda, sigmaY, Nparticles, x0, y, method = "mix")
# # runtime <- toc()
# # se <- (res_dac_mix$m[Time.step, ] - true_means[Time.step, ])^2
# # df <- data.frame(rbind(df, t(c(se, runtime))))
# # light
# tic()
# res_dac_light <- dac_time_lgssm_crossover(tau, lambda, sigmaY, Nparticles, x0, y, method = "light")
# runtime <- toc()
# se <- (res_dac_light$m[Time.step, ] - true_means[Time.step, ])^2
# df <- data.frame(rbind(df, t(c(se, runtime))))
# # adaptive light
# tic()
# res_dac_light_v2 <- dac_time_lgssm_crossover(tau, lambda, sigmaY, Nparticles, x0, y, method = "ada")
# runtime <- toc()
# se <- (res_dac_light_v2$m[Time.step, ] - true_means[Time.step, ])^2
# mean(se)
# df <- data.frame(rbind(df, t(c(se, runtime))))

# df$algo <- as.factor(c("dac", "mix", "light", "light_ada_temp", "dac", "mix", "light", "light_ada_temp"))
df$algo <- as.factor(c("dac", "light", "light_ada_temp"))
# df$mutation <- as.factor(c(0, 0, 0, 0, 1, 1, 1, 1))
df$mutation <- as.factor(c(0, 0, 0))

#df$algo <- as.factor(c("light_ada_temp", "light_ada_temp"))
#df$mutation <- as.factor(c(0, 1))
filename <- paste0("results/nocrossover_resampling_comparison_d", d, "N", Nparticles, "ID", ID)
write.csv(x=df, file=filename)
