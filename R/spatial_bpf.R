spatial_bpf <- function(sigmaX, y.error.var, y, Nparticles){
  d <- sqrt(nrow(y.error.var))
  x <- matrix(rnorm(Nparticles*d^2, sd = sqrt(sigmaX)), nrow = Nparticles, ncol = d^2)
  for (t in 2:(Time.step+1)){
    obs <- c(y[, , t-1])
    x <- x + matrix(rnorm(Nparticles*d^2, sd = sqrt(sigmaX)), nrow = Nparticles, ncol = d^2)
    lW <- -0.5*(nu+d^2)*log(1+apply(x, 1, FUN = function(z){ (obs - z)%*%y.error.var%*%(obs-z)})/nu)
    max.lW <- max(lW)
    W <- exp(lW - max(lW))
    W <- W/sum(W)
    ancestors <- stratified_resample(W, Nparticles)
    x <- x[ancestors, ]
  }
  return(x)
}
