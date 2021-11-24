mult_resample <- function(W, N){
  W <- c(W)
  # vector to store number of offsprings
  indices <- rep(0, times = N)

  # inverse transform
  un <- sort(runif(N))
  s <- W[1]
  m <- 1
  for(i in 1:N){
    while(s < un[i]){
      m <- m+1
      s <- s+W[m]
    }
    indices[i] <- m
  }
  return(indices)
}


stratified_resample <- function(W, N){
  W <- c(W)
  # vector to store number of offsprings
  indices <- rep(0, times = N)

  # inverse transform
  # partition (0, 1] into N disjoint intervals
  # and draw one uniform for every interval
  un <- (0:(N-1)+runif(N))/N
  s <- W[1]
  m <- 1
  for(i in 1:N){
    while(s < un[i]){
      m <- m+1
      s <- s+W[m]
    }
    indices[i] <- m
  }
  return(indices)
}
