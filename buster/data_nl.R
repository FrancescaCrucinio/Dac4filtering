# devtools::load_all("/storage/u1693998/Dac4filtering")

# ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
set.seed(1234*ID)

# dimension
d <- 4
# parameters
sigmaX <- 1
nu <- 10
delta <- 1
# precision of the observation equation
y.error.prec <- matrix(0, nrow = d^2, ncol = d^2)
diag(y.error.prec) <- 1
diag(y.error.prec[-1, ]) <- 1/4
vertical_neighbours <- ((0:d^2) * (d^2 + 1) + d+1)
y.error.prec[vertical_neighbours[(vertical_neighbours <= d^4)]] <- 1/4
y.error.prec[upper.tri(y.error.prec)] = t(y.error.prec)[upper.tri(y.error.prec)]

# number of time steps
Time.step <- 10

# get observations
nl_data <- nl_obs(d, sigmaX, nu, delta, y.error.prec, Time.step)

df_obs_iid <- data.frame(apply(nl_data$yiid,2,"c"))
df_obs_cov <- data.frame(apply(nl_data$y,2,"c"))
write.csv(x=df_obs_iid, file=paste0("/storage/u1693998/data/data_cov_nl_d", d, "ID", ID))
write.csv(x=df_obs_cov, file=paste0("/storage/u1693998/data/data_iid_nl_d", d, "ID", ID))
write.csv(x=data.frame(apply(nl_data$x,2,"c")), file=paste0("/storage/u1693998/data/data_truth_nl_d", d, "ID", ID))
