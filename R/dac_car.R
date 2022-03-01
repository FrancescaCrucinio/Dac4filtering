dac_car_lightweight <- function(history, obs, sigmaX, sigmaY, M = NULL){
  if(is.null(M)) {
    # number of samples for lightweight mixture (no adaptation)
    M <- ceiling(sqrt(Nparticles))
  }
  # dimension and number of particles
  d <- ncol(history[, , 1])
  Nparticles <- nrow(history[, , 1])
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # leaves
  # number of variables
  nv <- 1
  x <- matrix(0, nrow = Nparticles, ncol = d)
  lW <- matrix(0, nrow = Nparticles, ncol = d)
  W <- matrix(0, nrow = Nparticles, ncol = d)
  # history indices
  historyIndex <- array(1:Nparticles, dim = c(Nparticles, d, d))
  for (i in 1:nchild^nlevels){
    # propose
    x[, i] <- rowSums(history[, i:d, drop = FALSE, 1])/d + sqrt(sigmaX) * rnorm(Nparticles)
    # weights
    lW[, i] <- -0.5*(obs[i] - x[, i])^2/sigmaY - 0.5*log(2*pi*sigmaY)
    max.lW <- max(lW[, i])
    W[, i] <- exp(lW[, i] - max.lW)
    W[, i] <- W[, i]/sum(W[, i])
  }

  # loop over tree levels excluding leaves
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    nodes_dimension <- nchild^(u-1)

    # number of variables in each node
    nvNew <- nchild^u

    # updated particles
    xNew <- matrix(0, nrow = Nparticles, ncol = d)
    # updated history
    historyIndexNew <- array(0, dim = c(Nparticles, d, nodes))
    for (i in 1:nodes) {
      # get children indices
      ci <- child_indices(i, nvNew)
      if(u > 1){
        # mutation
        historyIndexNew[, , i] <- car_crossover(i, nodes, x, history, historyIndex, sigmaX)
      }
      else{ # at the leaf level all histories are the same
        historyIndexNew[, , i] <- historyIndex[, , i]
      }
      if(M == "adaptive"){
        # adaptive lightweight mixture resampling
        out <- car_adaptive_light(Nparticles, i, u, nv, nvNew, ci, lW, Nparticles, sigmaX, x, history[, , 1], historyIndex)
        # update after mixture resampling
        indices <- out$resampled_indices
        # update particles
        xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
        historyIndexNew[, , i] <- historyIndexNew[indices[, 1], , i]
        out$target_reached <- FALSE
        # if(!out$target_reached){
        #   # tempering
        #   tempering_out <- car_tempering(ci, i, nv, nvNew, sigmaX, sigmaY, obs, xNew, history[, , 1], historyIndex,
        #                                            historyIndexNew, out$resampled_particles_lW, Nparticles, 1 - 1e-05, 1/nodes_dimension)
        #   # update particles
        #   xNew <- tempering_out$x
        #   # update history
        #   historyIndexNew <- tempering_out$history_index_updated
        # }
        # else{
        #   # mcmc move
        #   updated_particles <- car_mcmc_move(x, ci, i, nv, nvNew, sigmaY, sigmaX, history[, , 1],
        #                                        historyIndex, obs, out$resampled_particles_lW, 1/nodes_dimension, 1)
        #   x <- updated_particles$x
        # }
      } else{
        out <- car_light(i, u, nv, nvNew, ci, W, Nparticles, M, sigmaX, x, history[, , 1], historyIndex)
        # update after mixture resampling
        indices <- out$resampled_indices
        xNew[, ci[1]:ci[2]] <- cbind(x[indices[, 1], ci[1]:(ci[1]+nv-1)], x[indices[, 2], (ci[1]+nv):ci[2]])
        historyIndexNew[, , i] <- historyIndexNew[indices[, 1], , i]
      }
    }
    x <- xNew
    nv <- nvNew
    historyIndex <- historyIndexNew
  }
  return(x)
}
