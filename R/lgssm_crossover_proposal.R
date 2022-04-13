crossover <- function(i, nodes, x, history, historyIndex, tau, lambda){
  # binary tree
  nchild <- 2
  Nparticles <- nrow(x)
  d <- ncol(x)
  historyIndexNew <- historyIndex[, , nchild*(i-1)+1]
  for (n in 1:Nparticles){
    crossover_point <- sample.int(d-1, 1)
    right_ancestor_coordinates <- cbind(historyIndex[n, (crossover_point+1):d, nchild*i], (crossover_point+1):d)
    left_ancestor_coordinates <- cbind(historyIndex[n, (crossover_point+1):d, nchild*(i-1)+1], (crossover_point+1):d)
    ft_ratio <- -0.5*(tau+lambda) * sum((x[n, (crossover_point+1):d] - (0.5*tau*history[right_ancestor_coordinates] + lambda*x[n, crossover_point:(d-1)])/(tau+lambda))^2 -
                  (x[n, (crossover_point+1):d] - (0.5*tau*history[left_ancestor_coordinates] + lambda*x[n, crossover_point:(d-1)])/(tau+lambda))^2)
    right_ancestor_coordinates <- cbind(historyIndex[n, crossover_point:(crossover_point+1), nchild*i], crossover_point:(crossover_point+1))
    left_ancestor_coordinates <- cbind(historyIndex[n, crossover_point:(crossover_point+1), nchild*(i-1)+1], crossover_point:(crossover_point+1))
    # turn into matrix to use as indices for history
    older_right_ancestor_coordinates <- matrix(c(historyIndex[n, crossover_point+1, nchild*i], crossover_point+1, 2), ncol = 3)
    older_left_ancestor_coordinates <- matrix(c(historyIndex[n, crossover_point+1, nchild*(i-1)+1], crossover_point+1, 2), ncol = 3)
    gamma_ratio <- -lambda * (history[right_ancestor_coordinates[1, , drop = FALSE]] - history[left_ancestor_coordinates[1, , drop = FALSE]]) *
                    (history[right_ancestor_coordinates[2, , drop = FALSE]] - history[left_ancestor_coordinates[2, , drop = FALSE]])
    # - 0.5*tau*history[older_right_ancestor_coordinates]/(tau+lambda) + 0.5*tau*history[older_left_ancestor_coordinates]/(tau+lambda))
    mh_ratio <- gamma_ratio + ft_ratio
    # accept/reject
    if(runif(1) <= exp(mh_ratio)){
      # accepted (only update the non-auxiliary history)
      historyIndexNew[n, ] <- c(historyIndex[n, 1:crossover_point, nchild*(i-1)+1], historyIndex[n, (crossover_point+1):d, nchild*i])
    }
  }
  return(historyIndexNew)
}
