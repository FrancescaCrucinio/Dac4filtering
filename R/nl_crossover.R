nl_crossover <- function(x, history, historyIndex_left, historyIndex_right, cir, cic, sigmaX, u,
                         covariance = FALSE, obs_old = NULL, tau = NULL){
  nu <- 10
  # crossover
  if(u > 1){
    if(covariance){
      merged_history <- nl_crossover_proposal_covariance(x, obs_old, history, historyIndex_left, historyIndex_right, cir, cic, sigmaX, nu, tau)
    }
    else {
      merged_history <- nl_crossover_proposal(x, history, historyIndex_left, historyIndex_right, cir, cic, sigmaX)
    }
  }
  else{ # at the leaf level all histories are the same
    merged_history <- historyIndex_left
  }
  return(merged_history)
}
