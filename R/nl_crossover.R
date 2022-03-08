nl_crossover <- function(x, history, historyIndex_left, historyIndex_right, nodes_row_left, nodes_row_right, nodes_col_left, nodes_col_right,
                         cir, cic, sigmaX, u){

  # crossover
  if(u > 1){
    merged_history <- nl_crossover_proposal(x, history, historyIndex_left, historyIndex_right, nodes_row_left, nodes_row_right,
                                            nodes_col_left, nodes_col_right, cir, cic, sigmaX)
  }
  else{ # at the leaf level all histories are the same
    merged_history <- historyIndex_left[, , ]
  }
  return(merged_history)
}
