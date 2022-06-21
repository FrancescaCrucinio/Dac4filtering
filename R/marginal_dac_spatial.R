marginal_dac_spatial <- function(history, obs, sigmaX, nu, tau, tau_diag, adaptive = FALSE){
  # dimension and number of particles
  d <- nrow(history)
  Nparticles <- dim(history)[3]
  y.error.var <- matrix(c(1, -0.25, -0.25, 0, -0.25, 1, 0, -0.25, -0.25, 0, 1, -0.25, 0, -0.25, -0.25, 1), nrow = 4)

  theta <- ceiling(sqrt(Nparticles))
  if(adaptive) theta <- NULL
  # tree topology
  nchild <- 2
  nlevels <- log2(d)
  # number of variables
  nv <- 1
  # leaves
  sample_from_past <- sample.int(Nparticles, Nparticles, replace = TRUE)
  x <- history[, , sample_from_past] + sqrt(sigmaX)*array(rnorm(Nparticles*d^2), dim = c(d, d, Nparticles))
  lW <- -0.5*(nu+1)*log(1+sweep(x, 1:2, obs)^2/nu)
  tmp <- matrix(x, ncol = 4, byrow = TRUE)
  lW2 <- -0.5*(nu+d^2)*log(1+apply(tmp, 1, FUN = function(z){ (c(obs) - z)%*%y.error.var%*%(c(obs)-z)})/nu)
  max.lW2 <- max(lW2)
  W2 <- exp(lW2 - max(lW2))
  W2 <- W2/sum(W2)
  ancestors <- stratified_resample(W2, Nparticles)
  print(paste(colMeans(tmp[ancestors, ])))
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
        out_top_merge <- marginal_spatial_merge(lW, x, obs, 2*i-1, 2*i-1, 2*j-1, 2*j,
                                           cir[, 1], cic[, 1], cir[, 1], cic[, 2],
                                           nv, nvNew, list("u" = u, "direction" = "h"), theta, tau, tau_diag, nu)
        x[cir[, 1], c(cic), ] <- out_top_merge$x
        print(paste(u, "h1"))
        # for (k in 1:d) {
        #   for (h in 1:d) {
        #     print(paste(mean(x[h, k, ])))
        #   }
        # }
        tmp <- matrix(x, ncol = 4, byrow = TRUE)
        lW2 <- -0.5*(nu+d^2)*log(1+apply(tmp, 1, FUN = function(z){ (c(obs) - z)%*%y.error.var%*%(c(obs)-z)})/nu)
        max.lW2 <- max(lW2)
        W2 <- exp(lW2 - max(lW2))
        W2 <- W2/sum(W2)
        ancestors <- stratified_resample(W2, Nparticles)
        print(paste(colMeans(tmp[ancestors, ])))
        ### Step 2
        out_bottom_merge <- marginal_spatial_merge(lW, x, obs, 2*i, 2*i, 2*j-1, 2*j,
                                                               cir[, 2], cic[, 1], cir[, 2], cic[, 2],
                                                               nv, nvNew, list("u" = u, "direction" = "h"), theta, tau, tau_diag, nu)
        x[cir[, 2], c(cic), ] <- out_bottom_merge$x
        print(paste(u, "h2"))
        # for (k in 1:d) {
        #   for (h in 1:d) {
        #     print(paste(mean(x[h, k, ])))
        #   }
        # }
        tmp <- matrix(x, ncol = 4, byrow = TRUE)
        lW2 <- -0.5*(nu+d^2)*log(1+apply(tmp, 1, FUN = function(z){ (c(obs) - z)%*%y.error.var%*%(c(obs)-z)})/nu)
        max.lW2 <- max(lW2)
        W2 <- exp(lW2 - max(lW2))
        W2 <- W2/sum(W2)
        ancestors <- stratified_resample(W2, Nparticles)
        print(paste(colMeans(tmp[ancestors, ])))
        #### VERTICAL MERGE ###
        out_merge <- marginal_spatial_merge(lW, x, obs, 2*i-1, 2*i, 2*j-1, 2*j-1,
                                                             cir[, 1], c(cic), cir[, 2], c(cic),
                                                             nv, nvNew, list("u" = u, "direction" = "v"), theta, tau, tau_diag, nu)
        x[c(cir), c(cic), ] <- out_merge$x
        print(paste(u, "v"))
        # for (k in 1:d) {
        #   for (h in 1:d) {
        #     print(paste(mean(x[h, k, ])))
        #   }
        # }
        tmp <- matrix(x, ncol = 4, byrow = TRUE)
        lW2 <- -0.5*(nu+d^2)*log(1+apply(tmp, 1, FUN = function(z){ (c(obs) - z)%*%y.error.var%*%(c(obs)-z)})/nu)
        max.lW2 <- max(lW2)
        W2 <- exp(lW2 - max(lW2))
        W2 <- W2/sum(W2)
        ancestors <- stratified_resample(W2, Nparticles)
        print(paste(colMeans(tmp[ancestors, ])))
      }
    }
    nv <- nvNew
  }
  return(x)
}
