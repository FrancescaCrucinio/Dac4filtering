nl_merge <- function(lW, obs, x, history, historyIndex, node_row_left, node_row_right,
                    node_col_left, node_col_right, cir_left, cic_left, cir_right, cic_right,
                    nv, nvNew, u_info, M, covariance = FALSE, tau = NULL){
  Nparticles <- dim(history)[3]
  d <- dim(history)[2]
  nchild <- 2
  sigmaX <- 1
  nu <- 10
  historyIndex_left <- historyIndex[, , , node_row_left, node_col_left]
  historyIndex_right <- historyIndex[, , , node_row_right, node_col_right]
  if(u_info$direction == "v"){
    u <- 2*u_info$u
  } else {
    u <- 2*u_info$u-1
  }
  if(u == 1){
    lW_left <- lW[node_row_left, node_col_left, ]
    lW_right <- lW[node_row_right, node_col_right, ]
  }
  if(covariance){
    if(is.null(M)){
      out <- nl_adaptive_light_covariance(Nparticles, u, obs, x, history, historyIndex_left, historyIndex_right, cir_left, cir_right, cic_left, cic_right,
                                          lW_left, lW_right, sigmaX, tau, nu, u_info)
      target_reached <- out$target_reached
    } else {
      target_reached <- TRUE
      out <- nl_light_covariance(u, obs, x, history, historyIndex_left, historyIndex_right, cir_left, cir_right, cic_left, cic_right,
                                 lW_left, lW_right, sigmaX, M, tau, nu)
    }
  } else{
    if(is.null(M)){
      out <- nl_adaptive_light(Nparticles, u, x, history, historyIndex_left, historyIndex_right, cir_right, cic_right,
                               lW_left, lW_right, sigmaX, u_info)
      target_reached <- out$target_reached
    } else {
      target_reached <- TRUE
      out <- nl_light(u, x, history, historyIndex_left, historyIndex_right, cir_right, cic_right,
                      lW_left, lW_right, sigmaX, M)
    }
  }

  indices <- out$resampled_indices
  if(u_info$direction == "h"){ # merge horizontally
    merged_x <- array(0, dim = c(nv, nvNew, Nparticles))
    merged_x[1:nv, 1:nv, ] <- x[cir_left, cic_left, indices[, 1], drop = FALSE]
    merged_x[1:nv, (nv+1):nvNew, ] <- x[cir_right, cic_right, indices[, 2], drop = FALSE]
  } else { # merge vertically
    merged_x <- array(0, dim = c(nvNew, nvNew, Nparticles))
    merged_x[1:nv, 1:nvNew, ] <- x[cir_left, cic_left, indices[, 1], drop = FALSE]
    merged_x[(nv+1):nvNew, 1:nvNew, ] <- x[cir_right, cic_right, indices[, 2], drop = FALSE]
  }
  return(list("x" = merged_x, "indices" = indices, "target_reached" = target_reached, "after_mix_lW" = out$resampled_particles_lW))
}
