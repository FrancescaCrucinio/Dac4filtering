# function returning index of first and last child of node i in level u
child_indices <- function(i, nv){
  rbind(nv*(i-1)+1, i*nv)
}
# get row and column indices for mixture resampling
rc_indices <- function(indices, dims){
  r <- indices -  trunc(indices/dims)*dims
  r[r == 0] <- dims
  c <- pmin(trunc((indices -1 )/dims)+1, dims)
  return(cbind(r, c))
}
# KS distance
ks_dist <- function(x, N){
  ks.test(x[1:N], x[N:length(x)])$statistic
}
# Wasserstein-1 distance
w1_dist <- function(x, N){
  wasserstein1d(x[1:N], x[N:length(x)])
}
