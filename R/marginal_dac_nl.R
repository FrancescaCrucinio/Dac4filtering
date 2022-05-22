marginal_dac_nl_lightweight <- function(history, obs, sigmaX, nu,
                               covariance = FALSE, obs_old = NULL, tau = NULL){
  # dimension and number of particles
  d <- nrow(history)
  Nparticles <- dim(history)[3]
  theta <- ceiling(sqrt(Nparticles))
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # number of variables
  nv <- 1
  x <- array(0, dim = c(d, d, Nparticles))
  lW <- array(0, dim = c(d, d, Nparticles))
  # leaves
  for (row in 1:d){
    for (col in 1:d){
      out_neighbours <- get_neighbours_weights(row, col, d)
      xMean <- sapply(1:Nparticles, marginal_sample_mixture, out_neighbours$mixture_weights,
                      out_neighbours$current_x_neighbours, history, simplify = TRUE)
      # weights
      x[row, col, ] <- xMean + sqrt(sigmaX*nlevels+1)*rnorm(Nparticles)
      lW[row, col, ] <- -0.5*(nu+1)*log(1+(x[row, col, ] - obs[row, col])^2/nu)
    }
  }
  # loop over tree levels excluding leaves
  for (u in 1:nlevels){
    # number of nodes at this level
    nodes <- nchild^(nlevels-u)
    # number of variables in each node
    nvNew <- nchild^u

    for (i in 1:nodes) {
      # row indices of children (column 1 = left child, column 2 = right child)
      cir <- matrix(c(((i-1)*nvNew+1):((2*i-1)*nv), ((2*i-1)*nv+1):(i*nvNew)), ncol = 2)
      for (j in 1:nodes) {
        # column indices of children (column 1 = left child, column 2 = right child)
        cic <- matrix(c(((j-1)*nvNew+1):((2*j-1)*nv), ((2*j-1)*nv+1):(j*nvNew)), ncol = 2)
        #### HORIZONTAL MERGE ###
        ### Step 1
        print("h1")
        out_top_merge <- marginal_nl_merge(lW, x, history, 2*i-1, 2*i-1, 2*j-1, 2*j, cir[, 1], cic[, 1],
                                  cir[, 1], cic[, 2], nv, nvNew, list("u" = u, "direction" = "h"), theta)
        x[cir[, 1], c(cic), ] <- out_top_merge$x
        ### Step 2
        print("h2")
        out_bottom_merge <- marginal_nl_merge(lW, x, history, 2*i, 2*i, 2*j-1, 2*j, cir[, 2], cic[, 1],
                                     cir[, 2], cic[, 2], nv, nvNew, list("u" = u, "direction" = "h"), theta)
        x[cir[, 2], c(cic), ] <- out_bottom_merge$x
        #### VERTICAL MERGE ###
        print("v")
        out_merge <- marginal_nl_merge(lW, x, history, 2*i-1, 2*i, 2*j-1, 2*j-1, cir[, 1], c(cic),
                              cir[, 2], c(cic), nv, nvNew, list("u" = u, "direction" = "v"), theta)
        x[c(cir), c(cic), ] <- out_merge$x
      }
    }
    nv <- nvNew
  }
  return(x)
}
