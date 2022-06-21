spatial_bpf <- function(history, sigmaX, nu, y.error.prec, obs, Nparticles){
  d <- sqrt(nrow(y.error.prec))
  x <- history + sqrt(sigmaX)*matrix(rnorm(Nparticles*d^2), nrow = Nparticles, ncol = d^2)
  lW <- -0.5*(nu+d^2)*log(1+apply(x, 1, FUN = function(z){ (obs - z)%*%y.error.prec%*%(obs-z)})/nu)
  max.lW <- max(lW)
  W <- exp(lW - max(lW))
  W <- W/sum(W)
  ancestors <- stratified_resample(W, Nparticles)
  x <- x[ancestors, ]
  return(x)
}
