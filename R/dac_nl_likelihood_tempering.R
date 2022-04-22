dac_nl_lightweight_likelihood_tempering <- function(history, obs, sigmaX, nu, M = NULL, covariance = FALSE, tempering = FALSE){
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
  W <- array(0, dim = c(d, d, Nparticles))
  # history indices
  historyIndex <- array(1:Nparticles, dim = c(Nparticles, d, d, d, d))
  # leaves
  beta <- (2-1)/(2-2^(-2*nlevels))
  for (col in 1:d){
    for (row in 1:d){
      out_neighbours <- get_neighbours_weights(row, col, d)
      xMean <- sapply(1:Nparticles, sample_mixture, out_neighbours$mixture_weights,
                      out_neighbours$current_x_neighbours, history, simplify = TRUE)
      # weights
      x[row, col, ] <- xMean + sqrt(sigmaX)*rnorm(Nparticles)
      lW[row, col, ] <- -beta*0.5*(nu+1)*log(1+(x[row, col, ] - obs[row, col])^2/nu)
      max.lW <- max(lW[row, col, ])
      W[row, col, ] <- exp(lW[row, col, ] - max.lW)
      W[row, col, ] <- W[row, col, ]/sum(W[row, col, ])
    }
  }
  # loop over tree levels excluding leaves
  for (u in 1:nlevels){

    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # nodes_dimension_h <- nchild^(u-1)
    # nodes_dimension_v <- nchild^(u)
    # number of variables in each node
    nvNew <- nchild^u
    # updated particles
    # xNew <- array(0, dim = c(d, d, Nparticles))
    for (i in 1:nodes) {
      # row indices of children (column 1 = left child, column 2 = right child)
      cir <- matrix(c(((i-1)*nvNew+1):((2*i-1)*nv), ((2*i-1)*nv+1):(i*nvNew)), ncol = 2)
      for (j in 1:nodes) {
        # column indices of children (column 1 = left child, column 2 = right child)
        cic <- matrix(c(((j-1)*nvNew+1):((2*j-1)*nv), ((2*j-1)*nv+1):(j*nvNew)), ncol = 2)
        #### HORIZONTAL MERGE ###
        ### Step 1
        # crossover
        historyIndexNew <- nl_crossover(x, history, historyIndex[, , , 2*i-1, 2*j-1], historyIndex[, , , 2*i-1, 2*j],
                                        cir[, 1], c(cic), sigmaX, u)
        # merge
        out_top_merge <- nl_merge_likelihood_tempering(lW, obs, x, history, historyIndex, 2*i-1, 2*i-1, 2*j-1, 2*j, cir[, 1], cic[, 1],
                                  cir[, 1], cic[, 2], nv, nvNew, list("u" = u, "direction" = "h"), M, covariance)
        x[cir[, 1], c(cic), ] <- out_top_merge$x
        historyIndex[, , , 2*i-1, 2*j-1] <- historyIndexNew[out_top_merge$indices[, 1], , , drop = FALSE]
        ### Step 2
        # crossover
        historyIndexNew <- nl_crossover(x, history, historyIndex[, , , 2*i, 2*j-1], historyIndex[, , , 2*i, 2*j],
                                           cir[, 2], c(cic), sigmaX, u)
        # merge
        out_bottom_merge <- nl_merge_likelihood_tempering(lW, obs, x, history, historyIndex, 2*i, 2*i, 2*j-1, 2*j, cir[, 2], cic[, 1],
                                     cir[, 2], cic[, 2], nv, nvNew, list("u" = u, "direction" = "h"), M, covariance)
        x[cir[, 2], c(cic), ] <- out_bottom_merge$x
        historyIndex[, , , 2*i, 2*j-1] <- historyIndexNew[out_bottom_merge$indices[, 1], , , drop = FALSE]
        #### VERTICAL MERGE ###
        # crossover
        historyIndexNew <- nl_crossover(x, history, historyIndex[, , , 2*i-1, 2*j-1], historyIndex[, , , 2*i, 2*j-1],
                                                    c(cir), c(cic), sigmaX, u+1)
        # merge
        out_merge <- nl_merge_likelihood_tempering(lW, obs, x, history, historyIndex, 2*i-1, 2*i, 2*j-1, 2*j-1, cir[, 1], c(cic),
                                                   cir[, 2], c(cic), nv, nvNew, list("u" = u, "direction" = "v"), M, covariance)
        x[c(cir), c(cic), ] <- out_merge$x
        historyIndexNew <- historyIndexNew[out_merge$indices[, 1], , ]
      }
    }
    # x <- xNew
    nv <- nvNew
  }
  return(x)
}
