library(ks)
node0 <- readRDS("test/node0.rds")
node1 <-  readRDS("test/node1.rds")
node2 <-  readRDS("test/node2.rds")
node3 <-  readRDS("test/node3.rds")
node4 <-  readRDS("test/node4.rds")
node5 <-  readRDS("test/node5.rds")
par(mfrow=c(4, 4))
for (i in 1:16) {
  kde0 <- kde(x = node0$x[, i], w = Nparticles*node0$W[, i])
  kde1 <- kde(x = node1[, i], eval.points = kde0$eval.points)
  kde2 <- kde(x = node2[, i], eval.points = kde0$eval.points)
  kde3 <- kde(x = node3[, i], eval.points = kde0$eval.points)
  kde4 <- kde(x = node4[, i], eval.points = kde0$eval.points)
  kde5 <- kde(x = node5[, i], eval.points = kde0$eval.points)
  truth <- dnorm(kde0$eval.points, mean = true_means[i], sd = sqrt(true_variances[i]))
  plot(kde0$eval.points, kde0$estimate, type = "l", col = "red", ylim = c(0, max(truth, kde0$estimate,
                                                                              kde1$estimate, kde2$estimate, kde3$estimate,
                                                                              kde4$estimate, kde5$estimate)))
                                                                                 # )))
  lines(kde0$eval.points, truth, type = "l", col = "black")
  lines(kde0$eval.points, kde1$estimate, type = "l", col = "blue")
  lines(kde0$eval.points, kde2$estimate, type = "l", col = "green")
  lines(kde0$eval.points, kde3$estimate, type = "l", col = "gray")
  lines(kde0$eval.points, kde4$estimate, type = "l", col = "brown")
  lines(kde0$eval.points, kde5$estimate, type = "l", col = "yellow")
}
par(mfrow=c(4, 4))

for (i in 17:32) {
  kde0 <- kde(x = node0$x[, i], w = Nparticles*node0$W[, i])
  kde1 <- kde(x = node1[, i], eval.points = kde0$eval.points)
  kde2 <- kde(x = node2[, i], eval.points = kde0$eval.points)
  kde3 <- kde(x = node3[, i], eval.points = kde0$eval.points)
  kde4 <- kde(x = node4[, i], eval.points = kde0$eval.points)
  kde5 <- kde(x = node5[, i], eval.points = kde0$eval.points)
  truth <- dnorm(kde0$eval.points, mean = true_means[i], sd = sqrt(true_variances[i]))
  plot(kde0$eval.points, kde0$estimate, type = "l", col = "red", ylim = c(0, max(truth, kde0$estimate,
                                                                                 kde1$estimate, kde2$estimate, kde3$estimate,
                                                                                 kde4$estimate, kde5$estimate)))
  # )))
  lines(kde0$eval.points, truth, type = "l", col = "black")
  lines(kde0$eval.points, kde1$estimate, type = "l", col = "blue")
  lines(kde0$eval.points, kde2$estimate, type = "l", col = "green")
  lines(kde0$eval.points, kde3$estimate, type = "l", col = "gray")
  lines(kde0$eval.points, kde4$estimate, type = "l", col = "brown")
  lines(kde0$eval.points, kde5$estimate, type = "l", col = "yellow")
}


