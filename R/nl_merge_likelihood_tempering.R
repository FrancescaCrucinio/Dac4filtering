nl_merge_likelihood_tempering <- function(lW, obs, x, history, historyIndex, node_row_left, node_row_right,
                     node_col_left, node_col_right, cir_left, cic_left, cir_right, cic_right, nv, nvNew, u_info, M, covariance = FALSE){
  Nparticles <- dim(history)[3]
  d <- dim(history)[2]
  nchild <- 2
  nlevels <- log2(d)
  sigmaX <- 1
  nu <- 10
  lW_left <- lW[node_row_left, node_col_left, ]
  lW_right <- lW[node_row_right, node_col_right, ]
  historyIndex_left <- historyIndex[, , , node_row_left, node_col_left]
  historyIndex_right <- historyIndex[, , , node_row_right, node_col_right]
  if(u_info$direction == "v"){
    u <- 2*u_info$u
  } else {
    u <- 2*u_info$u-1
  }
  if(covariance){
    if(is.null(M)){
      out <- nl_adaptive_light_covariance(Nparticles, u_info, obs, x, history, historyIndex_left, historyIndex_right, cir_left, cir_right, cic_left, cic_right,
                                          lW_left, lW_right, sigmaX, Nparticles, M, d, nchild^(u-1))
      target_reached <- out$target_reached
    } else {
      target_reached <- TRUE
      out <- nl_light_covariance(u, obs, x, history, historyIndex_left, historyIndex_right, cir_left, cir_right, cic_left, cic_right,
                                 lW_left, lW_right, sigmaX, Nparticles, M, d, nchild^(u-1))
    }
  } else{
    if(is.null(M)){
      out <- nl_adaptive_light(Nparticles, u_info, x, history, historyIndex_left, historyIndex_right, cir_right, cic_right,
                               lW_left, lW_right, sigmaX, M, d)
      target_reached <- out$target_reached
    } else {
      print(paste(u_info))
      beta_diff <- (-2^(-u)+2^(-u+1))/(2-2^(-2*nlevels))
      print(paste((2-2^(-u))/(2-2^(-2*nlevels))))
      target_reached <- TRUE
      out <- nl_light_likelihood_tempering(u, x, obs, history, historyIndex_left, historyIndex_right,
                                           cir_left, cir_right, cic_left, cic_right, lW_left, lW_right, sigmaX, nu, M, beta_diff)
    }
  }

  indices <- out$resampled_indices
  merged_x <- array(rbind(x[cir_left, cic_left, indices[, 1], drop = FALSE], x[cir_right, cic_right, indices[, 2], drop = FALSE]),
                    dim = c(nv, nvNew, Nparticles))
  return(list("x" = merged_x, "indices" = indices, "target_reached" = target_reached, "after_mix_lW" = out$resampled_particles_lW))
}
