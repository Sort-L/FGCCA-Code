GetMixCovMatrix = function(Ly.fun, 
                           Lt.fun, 
                           grid.fun,
                           Ly.mtv,
                           method="gam",
                           nbasis=10,
                           bw=NULL) {
  N <- length(Ly.fun)
  p.fun <- length(grid.fun)
  p.mtv <- dim(Ly.mtv)[2]
  cov.list <- list()
  # For each variable in multivariate block : smoothing
  for (j in seq_len(p.mtv)) {
    d.vec=rep(grid.fun, N)
    Lrcrcov <- lapply(seq_len(N), function(i){(Ly.fun[[i]] * Ly.mtv[i,j])})
    cov.list[[j]] <- GetMeanFun(Lrcrcov, Lt.fun, grid.fun, method, nbasis, bw)
  }
  CrCov <- do.call(rbind, cov.list)
  return(CrCov)
}