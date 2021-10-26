library(readr)
lgssm8_mse <- read_csv("buster/lgssm8_mse.csv")
lgssm8_runtime <- read_csv("buster/lgssm8_runtime.csv")

Time.step <- nrow(lgssm8_mse)
d <- ncol(lgssm8_mse)/2
# MSE
plot(1:Time.step, type = "l", rowMeans(lgssm8_mse[, 1:d]), col = "blue", xlab=" ", ylab=" ",
     ylim = c(0, max(rowMeans(lgssm8_mse))), cex = 2)
lines(1:Time.step, rowMeans(lgssm8_mse[, (d+1):(2*d)]), col = "red")
legend(1, 0.001, legend = c("dac", "dac-lw"), col=c("red", "blue"), lty=1, cex=1.5)
