crossover <- function(i, nodes, x, history, historyIndex, tau, lambda){
  # binary tree
  nchild <- 2
  Nparticles <- nrow(x)
  d <- ncol(x)
  # updated history
  historyIndexNew <- array(0, dim = c(Nparticles, d, nodes))
  for (n in 1:Nparticles){
    # sample crossover point
    crossover_point <- sample.int(d, 1)
    # accept/reject ratio
    mh_ratio <- -0.5*(tau+lambda) * sum( (x[n, (crossover_point+1):d]
                                          - 0.5*tau*history[historyIndex[n, (crossover_point+1):d, nchild*i], (crossover_point+1):d, 1]/(tau+lambda))^2 -
                                          (x[n, (crossover_point+1):d]
                                           - 0.5*tau*history[historyIndex[n, (crossover_point+1):d, nchild*(i-1)+1], (crossover_point+1):d, 1]/(tau+lambda))^2
                ) +
                0.5*tau*lambda * sum(x[n, crossover_point:(d-1)] *
                                       (history[historyIndex[n, (crossover_point+1):d, nchild*i], (crossover_point+1):d, 1]
                                        - history[historyIndex[n, (crossover_point+1):d, nchild*(i-1)+1], (crossover_point+1):d, 1])
                ) -
                lambda * (history[historyIndex[n, crossover_point, nchild*i], crossover_point, 1]
                          - history[historyIndex[n, crossover_point, nchild*(i-1)+1], crossover_point, 1]) *
                (history[historyIndex[n, crossover_point, nchild*i], crossover_point, 1]
                 - history[historyIndex[n, crossover_point, nchild*(i-1)+1], crossover_point, 1] -
                0.5*tau*history[historyIndex[n, crossover_point+1, nchild*i], crossover_point+1, 2] -
                  0.5*tau*history[historyIndex[n, crossover_point+1, nchild*(i-1)+1], crossover_point+1, 2])
    # accept/reject
    if(runif(1) <= exp(mh_ratio)){
      # accepted
      historyIndexNew[n, , i] <- c(historyIndex[n, 1:crossover_point, nchild*(i-1)+1], historyIndex[n, (crossover_point+1):d, nchild*i])
    }
  }
  return(historyIndexNew)
}
