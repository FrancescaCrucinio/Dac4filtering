# devtools::load_all("/storage/u1693998/Dac4filtering")
# ID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
ID <- 1
set.seed(1234*ID)
# get parameters
load("~/Documents/Dac4filtering/data/msv_parameters.RData")
# get data
y <- as.matrix(read.csv(file="data/synthetic_data_msv_y"))
true_x <- read.csv(file="data/synthetic_data_msv_x")
# dimension
p <- 26
r <- 4
d <- p + r

# number of time steps
Time.step <- 31
# number of particles
Nparticles <- 100

res_dac <- array(0, dim = c(Time.step+1, Nparticles, d))
res_dac[1, , ] <- mvrnorm(n = Nparticles, mu, diag(Sigma^2))
tic()
for (t in 1:Time.step) {
  res_dac[t+1, , ] <- marginal_dac_msv(res_dac[t, , ], y[t, ], mu, Phi, Lambda, Sigma, adaptive = FALSE)
  print(paste(t, t, t ,t, t, t))
}
toc()
#
# means.along <- function(a, i) {
#   n <- length(dim(a))
#   b <- aperm(a, c(seq_len(n)[-i], i))
#   rowMeans(b, dims = n - 1)
# }
# (means.along(res_dac, 2) - true_x[1:11, ])^2
# tmp = means.along(res_dac[1:123, ,], 2)
# colnames(tmp) <- x_levels
# library(tidyr)
# x_levels  <- vector(mode='character',length=30)
# for (i in 1:30) {
#   x_levels[i] <- paste0("X", i)
# }
# df <- gather(rbind(true_x[, ], tmp), "dim")
# df$time <- rep(1:123, times = 60)
# df$type <- rep(rep(c("true", "filter"), each = 123), times = 30)
#
# ggplot(data = df, aes(x = time, y = value, group = type, color= type)) +
#   geom_line(aes(group = type)) +
#   facet_wrap(~factor(dim, levels = x_levels), nrow = 15, ncol = 2, scales = "free", dir = "v") +
#   theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
#         axis.text = element_text(size=5),
#         # strip.text.x = element_blank(),
#         legend.title = element_blank(), legend.text=element_text(size=20),
#         text = element_text(size=5))


res <- fsvsample(y[1:31, ], factors = 4, draws = 10000, zeromean = TRUE, thin = 10, quiet = TRUE, restric = "auto")
cov_n <- rowMeans(covmat(res)[, , , 1], dims = 2)
library(reshape2)
melted_cormat <- melt(cov_n)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()
tmp <- colMeans(res_dac[31, , ])
estimated_cov <- rowMeans(res$facload, dims = 2) %*% diag(exp(tmp[(p+1):d])) %*% t(rowMeans(res$facload, dims = 2)) + diag(exp(tmp[1:p]))
