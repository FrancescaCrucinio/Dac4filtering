# devtools::load_all("/storage/u1693998/Dac4filtering")

ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
d <- 32
Time.step <- 100
Nparticles <- 1000
M <- 2*d
model = "lgssm"

output_stats(ID, d, Time.step, Nparticles, M, model)
