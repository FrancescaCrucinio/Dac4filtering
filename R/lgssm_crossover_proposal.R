crossover <- function(i, nodes, x, history, historyIndex, tau, lambda){
  # binary tree
  nchild <- 2
  Nparticles <- nrow(x)
  d <- ncol(x)
  historyIndexNew <- historyIndex[, , nchild*(i-1)+1]
  for (n in 1:Nparticles){
    crossover_point <- sample.int(d-1, 1)
    right_ancestor_coordinates <- cbind(historyIndex[n, (crossover_point+1):d, nchild*i], (crossover_point+1):d, rep(1, times = length((crossover_point+1):d)))
    left_ancestor_coordinates <- cbind(historyIndex[n, (crossover_point+1):d, nchild*(i-1)+1], (crossover_point+1):d, rep(1, times = length((crossover_point+1):d)))
    ft_ratio <- -0.5*(tau+lambda) * sum((x[n, (crossover_point+1):d] - 0.5*tau*history[right_ancestor_coordinates]/(tau+lambda))^2 -
                  (x[n, (crossover_point+1):d] - 0.5*tau*history[left_ancestor_coordinates]/(tau+lambda))^2) + 0.5*tau*lambda * sum(x[n, crossover_point:(d-1)] *
                  (history[right_ancestor_coordinates] - history[historyIndex[left_ancestor_coordinates], (crossover_point+1):d, 1]))
    older_right_ancestor_coordinates <- c(historyIndex[n, crossover_point+1, nchild*i], crossover_point+1, 2)
    older_left_ancestor_coodinates <- c(historyIndex[n, crossover_point+1, nchild*(i-1)+1], crossover_point+1, 2)
    gamma_ratio <- -lambda * (history[right_ancestor_coordinates] - history[left_ancestor_coordinates]) *
                    (history[right_ancestor_coordinates] - history[left_ancestor_coordinates] -
                    0.5*tau*history[older_right_ancestor_coordinates] - 0.5*tau*history[older_left_ancestor_coodinates])
    mh_ratio <- gamma_ratio + ft_ratio
    # accept/reject
    if(runif(1) <= exp(mh_ratio)){
      # accepted (only update the non-auxiliary history)
      historyIndexNew[n, ] <- c(historyIndex[n, 1:crossover_point, nchild*(i-1)+1], historyIndex[n, (crossover_point+1):d, nchild*i])
    }
  }
  return(historyIndexNew)
}
