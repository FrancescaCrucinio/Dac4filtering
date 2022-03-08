dac_nl_lightweight <- function(history, obs, sigmaX, nu, M = NULL){
  if(is.null(M)) {
    # number of samples for lightweight mixture (no adaptation)
    M <- ceiling(sqrt(Nparticles))
  }
  # dimension and number of particles
  d <- nrow(history[, , , 1])
  Nparticles <- dim(history[, , , 1])[3]
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # leaves
  # number of variables
  nv <- 1
  x <- array(0, dim = c(d, d, Nparticles))
  lW <- array(0, dim = c(d, d, Nparticles))
  # history indices
  # historyIndex <- array(1:Nparticles, dim = c(Nparticles, d, d))
  # leaves
  for (col in 1:d){
    for (row in 1:d){
      out_neighbours <- get_neighbours_weights(row, col, d)
      xMean <- sapply(1:Nparticles, sample_mixture, out_neighbours$mixture_weights,
                      out_neighbours$current_x_neighbours, history[, , , 1], simplify = TRUE)
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
    # historyIndexNew <- array(0, dim = c(Nparticles, d^2, nodes))
    for (i in 1:nodes) {
      # row indices of children (colum 1 = left child, column 2 = right child)
      cir <- matrix(c(((i-1)*nvNew+1):((2*i-1)*nv), ((2*i-1)*nv+1):(i*nvNew)), ncol = 2)
      for (j in 1:nodes) {
        # column indices of children (colum 1 = left child, column 2 = right child)
        cic <- matrix(c(((j-1)*nvNew+1):((2*j-1)*nv), ((2*j-1)*nv+1):(j*nvNew)), ncol = 2)
        #### HORIZONTAL MERGE ###
        ### Step 1
        left_child <- x[cir[, 1], cic[, 1], , drop = FALSE]
        right_child <- x[cir[, 1], cic[, 2], , drop = FALSE]
        lW_left <- lW[2*i-1, 2*j-1, ]
        lW_right <- lW[2*i-1, 2*j, ]
        out <- nl_light(u, left_child, right_child, lW_left, lW_right, Nparticles, M)
        indices <- out$resampled_indices
        xTop <- array(rbind(left_child[, ,indices[, 1], drop = FALSE], right_child[, , indices[, 2], drop = FALSE]), dim = c(nv, nvNew, Nparticles))
        ### Step 2
        left_child <- x[cir[, 2], cic[, 1], , drop = FALSE]
        right_child <- x[cir[, 2], cic[, 2], , drop = FALSE]
        lW_left <- lW[2*i, 2*j-1, ]
        lW_right <- lW[2*i, 2*j, ]
        out <- nl_light(u, left_child, right_child, lW_left, lW_right, Nparticles, M)
        indices <- out$resampled_indices
        xBottom <- array(rbind(left_child[, ,indices[, 1], drop = FALSE], right_child[, , indices[, 2], drop = FALSE]), dim = c(nv, nvNew, Nparticles))
        #### VERTICAL MERGE ###
        out <- nl_light(u+1, xTop, xBottom, lW_left, lW_right, Nparticles, M)
        indices <- out$resampled_indices
        xNew[c(cir), c(cic), ] <- array(rbind(xTop[, , indices[, 1], drop = FALSE], xBottom[, , indices[, 2], drop = FALSE]), dim = c(nvNew, nvNew, Nparticles))
      }
    }
    # historyIndex <- historyIndexNew
    x <- xNew
    nv <- nvNew
  }
  return(x)
}


