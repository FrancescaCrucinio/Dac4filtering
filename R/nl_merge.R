nl_merge <- function(lW, x, history, historyIndex, node_row_left, node_row_right,
                                node_col_left, node_col_right, cir_left, cic_left, cir_right, cic_right, nv, nvNew, u, M, covariance = FALSE){
  Nparticles <- dim(history)[3]
  d <- dim(history)[2]
  nchild <- 2
  lW_left <- lW[node_row_left, node_col_left, ]
  lW_right <- lW[node_row_right, node_col_right, ]
  historyIndex_left <- historyIndex[, , , node_row_left, node_col_left]
  historyIndex_right <- historyIndex[, , , node_row_right, node_col_right]
  if(covariance){
    if(is.null(M)){
      out <- nl_adaptive_light_covariance(Nparticles, u, x, history, historyIndex_left, historyIndex_right, cir_left, cir_right, cic_left, cic_right,
                               lW_left, lW_right, sigmaX, Nparticles, M, d, nchild^(u-1))
    } else {
      out <- nl_light_covariance(u, x, history, historyIndex_left, historyIndex_right, cir_left, cir_right, cic_left, cic_right,
                      lW_left, lW_right, sigmaX, Nparticles, M, d, nchild^(u-1))
    }
  } else{
    if(is.null(M)){
      out <- nl_adaptive_light(Nparticles, u, x, history, historyIndex_left, historyIndex_right, cir_right, cic_right,
                               lW_left, lW_right, sigmaX, Nparticles, M, d)
    } else {
      out <- nl_light(u, x, history, historyIndex_left, historyIndex_right, cir_right, cic_right,
                      lW_left, lW_right, sigmaX, Nparticles, M, d)
    }
  }

  indices <- out$resampled_indices
  merged_x <- array(rbind(x[cir_left, cic_left, indices[, 1], drop = FALSE], x[cir_right, cic_right, indices[, 2], drop = FALSE]),
                    dim = c(nv, nvNew, Nparticles))
  return(list("x" = merged_x, "indices" = indices))
}
