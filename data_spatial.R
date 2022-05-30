# Get data for spatial SSM
devtools::load_all("/storage/u1693998/Dac4filtering")
ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
set.seed(1234*ID)

# dimension
d <- 32
# parameters
sigmaX <- 1
nu <- 10
tau_diag <- 1
tau <- -0.25
# number of time steps
Time.step <- 10

# get observations
spatial_data <- spatial_obs(d, sigmaX, nu, tau, tau_diag, Time.step)

write.csv(x=data.frame(data.frame(apply(spatial_data$y,2,"c"))), file=paste0("/storage/u1693998/data/data_spatial_tau_", -tau, "d", d, "ID", ID))
write.csv(x=data.frame(apply(spatial_data$x,2,"c")), file=paste0("/storage/u1693998/data/data_truth_spatial_tau_", -tau, "d", d, "ID", ID))
