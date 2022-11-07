set.seed(1234)
# get true parameters
load("~/Documents/Dac4filtering/data/msv_parameters.RData")
# get data
y <- as.matrix(read.csv(file="data/synthetic_data_msv_y"))
# run factorstochvol on simulated data
res_factorstochvol <- fsvsample(y, factors = 4, draws = 10000, zeromean = TRUE, thin = 10, quiet = TRUE, restric = "auto")
corr_rep_factorstochvol <- array(apply(covmat(res_factorstochvol)[, , , 1], 3, cov2cor), c(26, 26, 1000))
cor_factorstochvol <- rowMeans(corr_rep_factorstochvol, dims = 2)
colnames(cor_factorstochvol) <- colnames(y)
rownames(cor_factorstochvol) <- colnames(y)
out <- corrplot(cor_factorstochvol, order = 'hclust')
# pdf(file = "synthetic_corrplot_factorstochvol.pdf")
# corrplot(cor_factorstochvol, order = 'hclust')
# dev.off()
Lambda_est <- rowMeans(res_factorstochvol$facload, dims = 2)
currency_order <- rownames(out$corr)

# dac results
Time.step <- 122
p <- 26
d <- 30
N <- "1000"
cor_dac_true_parameters <- array(0, dim = c(p, p, 50))
cor_dac_estimated_parameters <- array(0, dim = c(p, p, 50))
for (id in 1:50){
  filename <- paste0("data/msv/means_dac_msv_syn_d30N", N, "ID", id, sep = "")
  df <- read.csv(filename, row.names = 1)
  filtering_means <- df[Time.step, ]
  cor_dac_true_parameters[, , id] <- cov2cor(Lambda %*% diag(exp(filtering_means[(p+1):d])) %*% t(Lambda) + diag(exp(filtering_means[1:p])))
  cor_dac_estimated_parameters[, , id] <- cov2cor(Lambda_est %*% diag(exp(filtering_means[(p+1):d])) %*% t(Lambda_est) + diag(exp(filtering_means[1:p])))
}
cor_dac_true_parameters_mean <- rowMeans(cor_dac_true_parameters, dim = 2)
colnames(cor_dac_true_parameters_mean) <- colnames(y)
rownames(cor_dac_true_parameters_mean) <- colnames(y)
cor_dac_true_parameters_mean <- cor_dac_true_parameters_mean[currency_order, currency_order]
filename <- paste0("synthetic_corrplot_dac_true_parameters_N", N, ".pdf", sep = "")
out_dac <- corrplot(cor_dac_true_parameters_mean, order = 'original')
# pdf(file = filename)
# corrplot(cor_dac_true_parameters_mean, order = 'original')
# dev.off()

cor_dac_estimated_parameters_mean <- rowMeans(cor_dac_estimated_parameters, dim = 2)
colnames(cor_dac_estimated_parameters_mean) <- colnames(y)
rownames(cor_dac_estimated_parameters_mean) <- colnames(y)
cor_dac_estimated_parameters_mean <- cor_dac_estimated_parameters_mean[currency_order, currency_order]
out_dac2 <- corrplot(cor_dac_estimated_parameters_mean, order = 'original')
filename <- paste0("synthetic_corrplot_dac_estimated_parameters_N", N, ".pdf", sep = "")
# pdf(file = filename)
# out_dac2
# dev.off()


true_x <- as.matrix(read.csv(file="data/synthetic_data_msv_x"))[Time.step+1, ]
true_cov <- cov2cor(Lambda %*% diag(exp(true_x[(p+1):d])) %*% t(Lambda) + diag(exp(true_x[1:p])))
colnames(true_cov) <- colnames(y)
rownames(true_cov) <- colnames(y)
true_cov <- true_cov[currency_order, currency_order]
corrplot(true_cov, order = 'original')
