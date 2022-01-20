marginals <- matrix(0, nrow = 10^5, ncol = d)
for(i in 1:d){
  marginals[, i] <- rnorm(10^5, mean = true_means[Time.step, i], sd = sqrt(true_variances[Time.step, i]))
}


true_means <- data.matrix(data[data$type == "mean", 1:d])
true_variances <- data.matrix(data[data$type == "var", 1:d])
