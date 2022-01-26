car_crossover <- function(i, nodes, x, history, historyIndex, sigmaX){
  # binary tree
  nchild <- 2
  Nparticles <- nrow(x)
  d <- ncol(x)
  historyIndexNew <- historyIndex[, , nchild*(i-1)+1]
  for (n in 1:Nparticles){
    crossover_point <- sample.int(d-1, 1)
    left_ancestor_coordinates <- cbind(historyIndex[n, , nchild*(i-1)+1], 1:d, rep(1, times = d))
    right_ancestor_coordinates <- cbind(historyIndex[n, , nchild*i], 1:d, rep(1, times = d))
    crossed_over_history <- history[rbind(left_ancestor_coordinates[1:crossover_point, ], right_ancestor_coordinates[(crossover_point+1):d, ])]
    ft_ratio <- sum(-(x[n, ] - (c(0, cumsum(x[n, 1:(d-1)])) + rev(cumsum(rev(crossed_over_history))))/d)^2 +
                    (x[n, ] - (c(0, cumsum(x[n, 1:(d-1)])) + rev(cumsum(rev(history[left_ancestor_coordinates]))))/d)^2)/(2*sigmaX)

    gamma_ratio <- 0
    mh_ratio <- ft_ratio + gamma_ratio

    # accept/reject
    if(runif(1) <= exp(mh_ratio)){
      # accepted (only update the non-auxiliary history)
      historyIndexNew[n, ] <- c(historyIndex[n, 1:crossover_point, nchild*(i-1)+1], historyIndex[n, (crossover_point+1):d, nchild*i])
    }
  }
  return(historyIndexNew)
}
