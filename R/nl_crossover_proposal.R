nl_crossover_proposal <- function(x, history, historyIndex_left, historyIndex_right, cir, cic, sigmaX){
  d <- dim(history)[1]
  Nparticles <- dim(history)[3]
  # binary tree
  nchild <- 2

  merged_history <- historyIndex_left
  for (n in 1:Nparticles){
    crossover_point <- sample.int(d^2-1, 1)
    # crossover_point_col <- crossover_point -  trunc(crossover_point/d)*d
    # if(crossover_point_col == 0) crossover_point_col <- d
    # crossover_point_row <- pmin(trunc((crossover_point -1 )/d)+1, d)
    left_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_left[n, , ]))
    right_ancestor_coordinates <- cbind(1:d, rep(1:d, each = d), c(historyIndex_right[n, , ]))
    left_ancestor <- matrix(history[left_ancestor_coordinates], nrow = d)
    right_ancestor <- matrix(history[right_ancestor_coordinates], nrow = d)
    crossedover_ancestor <- matrix(c(left_ancestor[1:crossover_point], right_ancestor[(crossover_point+1):d^2]), nrow = d)
    ft_ratio <- 0
    for (col in cic) {
      for (row in cir) {
        out_neighbours <- get_neighbours_weights(row, col, d)
        valid_weights <- out_neighbours$mixture_weights[out_neighbours$mixture_weights>0]
        valid_current_neighbours <- out_neighbours$current_x_neighbours[out_neighbours$mixture_weights>0, ]
        ft_ratio <- ft_ratio + log(sum(valid_weights * dnorm(x[row, col, n], mean = crossedover_ancestor[valid_current_neighbours], sd = sqrt(sigmaX)))) -
                log(sum(valid_weights * dnorm(x[row, col, n], mean = left_ancestor[valid_current_neighbours], sd = sqrt(sigmaX))))
      }
    }
    # accept/reject
    if(runif(1) <= exp(ft_ratio)){
      # accepted (only update the non-auxiliary history)
      tmp_left <- historyIndex_left[n, , ]
      tmp_right <- historyIndex_right[n, , ]
      merged_history[n, , ] <- matrix(c(tmp_left[1:crossover_point], tmp_right[(crossover_point+1):d^2]),
                                      ncol = dim(historyIndex_left)[2])
    }
  }
  return(merged_history)
}
