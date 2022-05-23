marginal_nl_merge <- function(lW, x, history, node_row_left, node_row_right,
                     node_col_left, node_col_right, cir_left, cic_left, cir_right, cic_right,
                     nv, nvNew, u_info, theta){
  Nparticles <- dim(history)[3]
  d <- dim(history)[2]
  sigmaX <- 1
  if(u_info$u == 1 & u_info$direction == "h"){
    lW_left <- lW[node_row_left, node_col_left, ]
    lW_right <- lW[node_row_right, node_col_right, ]
  }
  if(is.null(theta)){
    out <- marginal_nl_light_fast_adaptive(Nparticles, u_info, x, history, cir_left, cir_right, cic_left, cic_right,
                                  lW_left, lW_right, sigmaX)
  } else {
    out <- marginal_nl_light_fast(u_info, x, history, cir_left, cir_right, cic_left, cic_right,
                                  lW_left, lW_right, sigmaX, theta)
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
  return(list("x" = merged_x, "indices" = indices, "after_mix_lW" = out$resampled_particles_lW))
}
