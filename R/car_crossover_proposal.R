car_crossover <- function(i, nodes, x, history, historyIndex, sigmaX){
  # binary tree
  nchild <- 2
  Nparticles <- nrow(x)
  d <- ncol(x)
  historyIndexNew <- historyIndex[, , nchild*(i-1)+1]
  for (n in 1:Nparticles){
    # sample crossover point
    crossover_point <- sample.int(d-1, 1)
    # accept/reject ratio
    # mh_ratio <- -((x[n, ] - (c(0, cumsum(x[n, ])[-d]) + )/d)^2 +
    #                 (x[n, ] - (c(0, cumsum(x[n, ])[-d]) + )/d)^2)/(2*sigmaX)
    #
    #
    # for (h in 1:d){
    #   mh_ratio <- mh_ratio -
    #
    #
    #     ((x[n, h] -d^{-1}(sum(x[n, 1:(h-1)]) + ifelse(h<=crossover_point, sum(history[historyIndex[n, h:crossover_point, nchild*(i-1)+1], h:crossover_point, 1], history[historyIndex[n, (crossover_point+1):d, nchild*i], (crossover_point+1):d, 1]),
    #                                                                        sum(history[historyIndex[n, h:d, nchild*i], h:d, 1]))))^2 +
    #     (x[n, h] -d^{-1}(sum(x[n, 1:(h-1)] + sum(history[historyIndex[n, h:d, nchild*(i-1)+1], h:d, 1]))))^2
    #     )/(2*sigmaX)
    # }

    # accept/reject
    if(runif(1) <= exp(mh_ratio)){
      # accepted (only update the non-auxiliary history)
      historyIndexNew[n, ] <- c(historyIndex[n, 1:crossover_point, nchild*(i-1)+1], historyIndex[n, (crossover_point+1):d, nchild*i])
    }
  }
  return(historyIndexNew)
}
