marginal_spatial_merge <- function(lW, x, obs, node_row_left, node_row_right,
                              node_col_left, node_col_right, cir_left, cic_left, cir_right, cic_right,
                              nv, nvNew, u_info, theta, tau, nu){
  Nparticles <- dim(x)[3]
  d <- dim(x)[2]
  if(u_info$u == 1 & u_info$direction == "h"){
    lW_left <- lW[node_row_left, node_col_left, ]
    lW_right <- lW[node_row_right, node_col_right, ]
  }
  if(is.null(theta)){
    out <- marginal_spatial_light_adaptive(Nparticles, u_info, x, obs, cir_left, cir_right, cic_left, cic_right,
                                           lW_left, lW_right, sigmaX, tau, nu)
  } else {
    out <- marginal_spatial_light(u_info, x, obs, cir_left, cir_right, cic_left, cic_right,
                                  lW_left, lW_right, sigmaX, theta, tau, nu)
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
