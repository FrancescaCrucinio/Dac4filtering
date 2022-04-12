dac_nl_lightweight <- function(history, obs, sigmaX, nu, M = NULL,
                               covariance = FALSE, tempering = FALSE, obs_old = NULL, tau = NULL){
  # dimension and number of particles
  d <- nrow(history)
  Nparticles <- dim(history)[3]
  if(is.null(M)) {
    # number of samples for lightweight mixture (no adaptation)
    M <- ceiling(sqrt(Nparticles))
  }
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # leaves
  # number of variables
  nv <- 1
  x <- array(0, dim = c(d, d, Nparticles))
  lW <- array(0, dim = c(d, d, Nparticles))
  # history indices
  historyIndex <- array(1:Nparticles, dim = c(Nparticles, d, d, d, d))
  # leaves
  for (col in 1:d){
    for (row in 1:d){
      out_neighbours <- get_neighbours_weights(row, col, d)
      xMean <- sapply(1:Nparticles, sample_mixture, out_neighbours$mixture_weights,
                      out_neighbours$current_x_neighbours, history, simplify = TRUE)
      # weights
      x[row, col, ] <- xMean + sqrt(sigmaX)*rnorm(Nparticles)
      lW[row, col, ] <- -0.5*(nu+1)*log(1+(x[row, col, ] - obs[row, col])^2/nu)
    }
  }
  # loop over tree levels excluding leaves
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    nodes_dimension_h <- nchild^(u-1)
    nodes_dimension_v <- nchild^(u)
    # number of variables in each node
    nvNew <- nchild^u

    # updated particles
    xNew <- array(0, dim = c(d, d, Nparticles))
    # updated history
    historyIndexNew <- array(0, dim = c(Nparticles, d, d, nodes, nodes))
    for (i in 1:nodes) {
      # row indices of children (column 1 = left child, column 2 = right child)
      cir <- matrix(c(((i-1)*nvNew+1):((2*i-1)*nv), ((2*i-1)*nv+1):(i*nvNew)), ncol = 2)
      for (j in 1:nodes) {
        # column indices of children (column 1 = left child, column 2 = right child)
        cic <- matrix(c(((j-1)*nvNew+1):((2*j-1)*nv), ((2*j-1)*nv+1):(j*nvNew)), ncol = 2)
        #### HORIZONTAL MERGE ###
        ### Step 1
        # crossover
        historyIndexTop <- nl_crossover(x, history, historyIndex[, , , 2*i-1, 2*j-1], historyIndex[, , , 2*i-1, 2*j],
                                        cir[, 1], c(cic), sigmaX, u, covariance, obs_old, tau)
        # merge
        out_top_merge <- nl_merge(lW, obs, x, history, historyIndex, 2*i-1, 2*i-1, 2*j-1, 2*j, cir[, 1], cic[, 1],
                                  cir[, 1], cic[, 2], nv, nvNew, list("u" = u, "direction" = "h"), M, covariance, tau)
        xNew[cir[, 1], c(cic), ] <- out_top_merge$x
        historyIndex[, , , 2*i-1, 2*j-1] <- historyIndexTop[out_top_merge$indices[, 1], , , drop = FALSE]
        if(tempering){
          lW_left <- lW[2*i-1, 2*j-1, ]
          lW_right <- lW[2*i-1, 2*j, ]
          # tempering
          tempering_out_top <- nl_tempering(Nparticles, 1-1e-05, x, history, historyIndex[, , , 2*i-1, 2*j-1], historyIndex[, , , 2*i-1, 2*j],
                                            obs, u, nu, out_top_merge$after_mix_lW, lW_left, lW_right, c(cir), c(cic),
                                            cir[, 1], cir[, 1], cic[, 1], cic[, 2], sqrt(1/nodes_dimension_h))
          # update particles
          xNew <- tempering_out_top$x
          # update history
          historyIndex[, , , 2*i-1, 2*j-1] <- tempering_out_top$history_index_updated
        }
        ### Step 2
        # crossover
        historyIndexBottom <- nl_crossover(x, history, historyIndex[, , , 2*i, 2*j-1], historyIndex[, , , 2*i, 2*j],
                                           cir[, 2], c(cic), sigmaX, u, covariance, obs_old, tau)
        # merge
        out_bottom_merge <- nl_merge(lW, obs, x, history, historyIndex, 2*i, 2*i, 2*j-1, 2*j, cir[, 2], cic[, 1],
                                     cir[, 2], cic[, 2], nv, nvNew, list("u" = u, "direction" = "h"), M, covariance, tau)
        xNew[cir[, 2], c(cic), ] <- out_bottom_merge$x
        historyIndex[, , , 2*i, 2*j-1] <- historyIndexBottom[out_bottom_merge$indices[, 1], , , drop = FALSE]
        if(tempering){
          lW_left <- lW[2*i, 2*j-1, ]
          lW_right <- lW[2*i, 2*j, ]
          # tempering
          tempering_out_bottom <- nl_tempering(Nparticles, 1-1e-05, x, history, historyIndex[, , , 2*i, 2*j-1], historyIndex[, , , 2*i, 2*j],
                                            obs, u, nu, out_bottom_merge$after_mix_lW, lW_left, lW_right, c(cir), c(cic),
                                            cir[, 2], cir[, 2], cic[, 1], cic[, 2], sqrt(1/nodes_dimension_h))
          # update particles
          xNew <- tempering_out_bottom$x
          # update history
          historyIndex[, , , 2*i, 2*j-1] <- tempering_out_bottom$history_index_updated
        }
        #### VERTICAL MERGE ###
        # crossover
        historyIndexNew[, , , i, j] <- nl_crossover(xNew, history, historyIndex[, , , 2*i-1, 2*j-1], historyIndex[, , , 2*i, 2*j-1],
                                                    c(cir), c(cic), sigmaX, u+1, covariance, obs_old, tau)
        # merge
        out_merge <- nl_merge(lW, obs, xNew, history, historyIndex, 2*i-1, 2*i, 2*j-1, 2*j-1, cir[, 1], c(cic),
                              cir[, 2], c(cic), nvNew, nvNew, list("u" = u, "direction" = "v"), M, covariance, tau)
        xNew[c(cir), c(cic), ] <- out_merge$x
        historyIndexNew[, , , i, j] <- historyIndexNew[out_merge$indices[, 1], , , i, j]
        if(tempering){
          # tempering
          tempering_out_merge <- nl_tempering(Nparticles, 1-1e-05, xNew, history, historyIndex[, , , 2*i-1, 2*j-1], historyIndex[, , , 2*i, 2*j-1],
                                            obs, u+1, nu, out_merge$after_mix_lW, lW_left, lW_right, c(cir), c(cic),
                                            cir[, 2], cir[, 1], c(cic), c(cic), sqrt(1/nodes_dimension_v))
          # update particles
          xNew <- tempering_out_merge$x
          # update history
          historyIndexNew[, , , i, j] <- tempering_out_merge$history_index_updated
        }
      }
    }
    historyIndex <- historyIndexNew
    x <- xNew
    nv <- nvNew
  }
  return(x)
}


dac_nl_adaptive_lightweight <- function(history, obs, sigmaX, nu, M = NULL, covariance = FALSE){
  # dimension and number of particles
  d <- nrow(history)
  Nparticles <- dim(history)[3]
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # leaves
  # number of variables
  nv <- 1
  x <- array(0, dim = c(d, d, Nparticles))
  lW <- array(0, dim = c(d, d, Nparticles))
  # history indices
  historyIndex <- array(1:Nparticles, dim = c(Nparticles, d, d, d, d))
  # leaves
  for (col in 1:d){
    for (row in 1:d){
      out_neighbours <- get_neighbours_weights(row, col, d)
      xMean <- sapply(1:Nparticles, sample_mixture, out_neighbours$mixture_weights,
                      out_neighbours$current_x_neighbours, history, simplify = TRUE)
      # weights
      x[row, col, ] <- xMean + sqrt(sigmaX)*rnorm(Nparticles)
      lW[row, col, ] <- -0.5*(nu+1)*log(1+(x[row, col, ] - obs[row, col])^2/nu)
    }
  }
  # loop over tree levels excluding leaves
  for (u in 1:nlevels){

    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u

    # updated particles
    xNew <- array(0, dim = c(d, d, Nparticles))
    # updated history
    historyIndexNew <- array(0, dim = c(Nparticles, d, d, nodes, nodes))
    for (i in 1:nodes) {
      # row indices of children (colum 1 = left child, column 2 = right child)
      cir <- matrix(c(((i-1)*nvNew+1):((2*i-1)*nv), ((2*i-1)*nv+1):(i*nvNew)), ncol = 2)
      for (j in 1:nodes) {
        # column indices of children (colum 1 = left child, column 2 = right child)
        cic <- matrix(c(((j-1)*nvNew+1):((2*j-1)*nv), ((2*j-1)*nv+1):(j*nvNew)), ncol = 2)
        #### HORIZONTAL MERGE ###
        ### Step 1
        # crossover
        historyIndexTop <- nl_crossover(x, history, historyIndex[, , , 2*i-1, 2*j-1], historyIndex[, , , 2*i-1, 2*j],
                                        cir[, 1], c(cic), sigmaX, u)
        # merge
        out_top_merge <- nl_merge(lW, obs, x, history, historyIndex, 2*i-1, 2*i-1, 2*j-1, 2*j, cir[, 1], cic[, 1],
                                  cir[, 1], cic[, 2], nv, nvNew, list("u" = u, "direction" = "h"), M, covariance)
        xNew[cir[, 1], c(cic), ] <- out_top_merge$x
        historyIndex[, , , 2*i-1, 2*j-1] <- historyIndexTop[out_top_merge$indices[, 1], , , drop = FALSE]
        ### Step 2
        # crossover
        historyIndexBottom <- nl_crossover(x, history, historyIndex[, , , 2*i, 2*j-1], historyIndex[, , , 2*i, 2*j],
                                           cir[, 2], c(cic), sigmaX, u)
        # merge
        out_bottom_merge <- nl_merge(lW, obs, x, history, historyIndex, 2*i, 2*i, 2*j-1, 2*j, cir[, 2], cic[, 1],
                                     cir[, 2], cic[, 2], nv, nvNew, list("u" = u, "direction" = "h"), M, covariance)
        xNew[cir[, 2], c(cic), ] <- out_bottom_merge$x
        historyIndex[, , , 2*i, 2*j-1] <- historyIndexBottom[out_bottom_merge$indices[, 1], , , drop = FALSE]
        #### VERTICAL MERGE ###
        # crossover
        historyIndexNew[, , , i, j] <- nl_crossover(xNew, history, historyIndex[, , , 2*i-1, 2*j-1], historyIndex[, , , 2*i, 2*j-1],
                                                    c(cir), c(cic), sigmaX, u+1)
        # merge
        out_merge <- nl_merge(lW, obs, xNew, history, historyIndex, 2*i-1, 2*i, 2*j-1, 2*j-1, cir[, 1], c(cic),
                              cir[, 2], c(cic), nvNew, nvNew, list("u" = u, "direction" = "v"), M, covariance)
        xNew[c(cir), c(cic), ] <- out_merge$x
        historyIndexNew[, , , i, j] <- historyIndexNew[out_merge$indices[, 1], , , i, j]
      }
    }
    historyIndex <- historyIndexNew
    x <- xNew
    nv <- nvNew
  }
  return(x)
}


