library("factorstochvol")
library("zoo")
# data
load("/Users/francescacrucinio/Documents/Dac4filtering/data/exrates.RData")
dat$date <- as.POSIXct(rownames(dat), format = "%Y-%m-%d")
# get data unti January 2013
exrates_to2013 <- dat[dat$date <= "2013-01-04", ]
first_semester2013 <- dat[dat$date > "2013-01-04" & dat$date < "2013-07-01", ]
m <- 26 # currencies
n <- dim(exrates_to2013)[1] # days
y <- 100 * logret(tail(exrates_to2013[, seq_len(m)], n + 1))
plot(zoo(y, order.by = tail(exrates_to2013$date, n)), main = "", xlab = "Time")
# run MCMC
set.seed(1234)
res <- fsvsample(y, factors = 4, draws = 10000, zeromean = TRUE, thin = 10, quiet = TRUE, restric = "auto")
# diagnostics
dim(cov_n <- covmat(res))
logdet <- function (x) log(det(x))
logdet_n <- apply(cov_n[,,,1], 3, logdet)
ts.plot(logdet_n)
acf(logdet_n, main = "")
# # beta parameter
# rowMeans(res$beta)
# mu parameter
mu <- rowMeans(res$para[1, , ])
# Phi parameter
Phi <- rowMeans(res$para[2, , ])
# sigma parameter
Sigma <- rowMeans(res$para[3, , ])
# lambda parameter
Lambda <- rowMeans(res$facload, dims = 2)
save(mu, Phi, Sigma, Lambda, list = c("mu", "Phi", "Sigma", "Lambda"), file = "msv_parameters.RData")

