nl_merge <- function(lW, xleft, xright, historyIndex, node_row_left, node_row_right,
                                node_col_left, node_col_right, nv, nvNew, u, M, Nparticles){
  # merge
  lW_left <- lW[node_row_left, node_col_left, ]
  lW_right <- lW[node_row_right, node_col_right, ]
  out <- nl_light(u, xleft, xright, lW_left, lW_right, Nparticles, M)
  indices <- out$resampled_indices
  merged_x <- array(rbind(xleft[, ,indices[, 1], drop = FALSE], xright[, , indices[, 2], drop = FALSE]),
                    dim = c(nv, nvNew, Nparticles))
  return(list("x" = merged_x, "history" = historyIndex))
}
