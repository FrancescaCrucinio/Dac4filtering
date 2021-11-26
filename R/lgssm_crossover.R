crossover <- function(crossover_point, ci, z, history, tau, lambda){
  # merge the x on the two children
  # merged_x <- c(x[n1, ci[1]:(ci[1]+nv-1)], x[n2, (ci[1]+nv)]
  # accept/reject
  mh_ratio <- -0.5*(tau+lambda) * sum((z[(crossover_point+1):d] - 0.5*tau*history[n2, (crossover_point+1):d, 1]/(tau+lambda))^2 -
                                      (z[(crossover_point+1):d] - 0.5*tau*history[n1, (crossover_point+1):d, 1]/(tau+lambda))^2) -
              -0.5*tau*lambda * z[(crossover_point):d] * (history[n2, (crossover_point+1):d, 1] - history[n1, (crossover_point+1):d, 1])
              -lambda * (history[n2, crossover_point, 1] - history[n1, crossover_point, 1]) *
                (history[n2, crossover_point, 1] - history[n1, crossover_point, 1] -
                   0.5*tau*history[n2, crossover_point+1, 2] - 0.5*tau*history[n1, crossover_point+1, 2])


  if(runif(1) <= exp(mh_ratio)){
    history[n1, (crossover_point+1):ci[2], 1] <- history[n2, (crossover_point+1):ci[2], 1]
  }
  return(history)
}
