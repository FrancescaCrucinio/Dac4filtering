# data
load("/Users/francescacrucinio/Documents/Dac4filtering/data/exrates.RData")
dat$date <- as.POSIXct(rownames(dat), format = "%Y-%m-%d")
### data until January 2013
exrates_to2013 <- dat[dat$date <= "2013-01-04", ]
m <- 26 # currencies
n <- dim(exrates_to2013)[1] # days
y <- 100 * logret(tail(exrates_to2013[, seq_len(m)], n + 1))
# run MCMC
set.seed(1234)
res <- fsvsample(y, factors = 4, draws = 10000, zeromean = TRUE, thin = 10, quiet = TRUE, restric = "auto")
# diagnostics
dim(cov_n <- covmat(res))
logdet <- function (x) log(det(x))
logdet_n <- apply(cov_n[,,,1], 3, logdet)
ts.plot(logdet_n)
acf(logdet_n, main = "")
# mu parameter
mu <- rowMeans(res$para[1, , ])
# Phi parameter
Phi <- rowMeans(res$para[2, , ])
# sigma parameter
Sigma <- rowMeans(res$para[3, , ])
# lambda parameter
Lambda <- rowMeans(res$facload, dims = 2)
save(mu, Phi, Sigma, Lambda, list = c("mu", "Phi", "Sigma", "Lambda"), file = "msv_parameters.RData")

### first semester of 2013
first_semester2013 <- dat[dat$date > "2013-01-04" & dat$date <= "2013-07-01", ]
write.csv(x=first_semester2013[, 1:26], file="data/real_data_msv_y", row.names = FALSE)
# run factorstochvol on last semester
n <- dim(first_semester2013)[1] # days
y_semester <- 100 * logret(tail(first_semester2013[, seq_len(m)], n + 1))
res_semester <- fsvsample(y, factors = 4, draws = 10000, zeromean = TRUE, thin = 10, quiet = TRUE, restric = "auto")
cov_n_semester <- rowMeans(covmat(res_semester)[, , , 1], dims = 2)
colnames(cov_n_semester) <- colnames(dat)[1:26]
rownames(cov_n_semester) <- colnames(dat)[1:26]
write.csv(x=data.frame(cov_n_semester), file="data/real_data_msv_cov_mat", row.names = FALSE)

