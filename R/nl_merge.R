nl_merge <- function(lW, x, history, historyIndex, node_row_left, node_row_right,
                                node_col_left, node_col_right, cir_left, cic_left, cir_right, cic_right, nv, nvNew, u, M){
  Nparticles <- dim(history)[3]
  d <- dim(history)[1]
  lW_left <- lW[node_row_left, node_col_left, ]
  lW_right <- lW[node_row_right, node_col_right, ]
  historyIndex_left <- historyIndex[, , , node_row_left, node_col_left]
  historyIndex_right <- historyIndex[, , , node_row_right, node_col_right]

  out <- nl_light(u, x, history, historyIndex_left, historyIndex_right, cir_right, cic_right,
                  lW_left, lW_right, sigmaX, Nparticles, M, d)
  indices <- out$resampled_indices
  merged_x <- array(rbind(x[cir_left, cic_left, indices[, 1], drop = FALSE], x[cir_right, cic_right, indices[, 2], drop = FALSE]),
                    dim = c(nv, nvNew, Nparticles))
  return(list("x" = merged_x, "history" = historyIndex))
}
