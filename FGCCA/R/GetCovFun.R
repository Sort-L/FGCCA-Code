GetCovFun <- function(Ly, 
                      Lt, 
                      grid,
                      smoothing="gam",
                      nbasis=10,
                      bw=NULL) {
  obsGrid <- sort(unique(unlist(Lt)))
  N <- length(Ly)
  p <- length(obsGrid)
  cov.sum = cov.count = cov.mean = matrix(0, p, p)
  for (i in seq_len(N)) {
    nas.points = which(!is.na(Ly[[i]]))
    obs.points = which(obsGrid %in% Lt[[i]][nas.points])
    cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] + 1
    cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] + 
      tcrossprod(Ly[[i]][nas.points])
  }
  cov.mean = ifelse(cov.count == 0, NA, cov.sum/cov.count)
  diag.cov = diag(cov.mean)
  diag.cnt = diag(cov.count)
  diag(cov.mean) = NA
  
  # Removing zeros
  row.vec = rep(obsGrid, each=p)
  col.vec = rep(obsGrid, p)
  cov.vec = as.vector(cov.mean)
  cnt.vec = as.vector(matrix(1, p, p))
  
  row.vec.pred = rep(grid, each=length(grid))
  col.vec.pred = rep(grid, length(grid))
  
  if (smoothing == "gam") {
    CrCov = matrix(predict(mgcv::gam(cov.vec ~ te(row.vec, col.vec, k=nbasis), weights = cnt.vec),
                           newdata = data.frame(row.vec = row.vec.pred, col.vec = col.vec.pred)), nrow=length(grid), ncol=length(grid))
    DiagCov = as.vector(predict(mgcv::gam(diag.cov ~ s(obsGrid, k=nbasis), weights = diag.cnt),
                                newdata = data.frame(obsGrid = grid)))
    CrCov = (CrCov + t(CrCov)) / 2
  } else {
    nonzeros.vals = which(!is.na(cov.vec))
    pairs.df = cbind(row.vec[nonzeros.vals], col.vec[nonzeros.vals])
    CrCov = fdapace::Lwls2D(c(bw, bw), kern="gauss", xin=pairs.df, yin=cov.vec[nonzeros.vals],
                   win = cnt.vec[nonzeros.vals], xout1=grid, xout2=grid)
    DiagCov = fdapace::Lwls1D(bw, kern="gauss", xin=obsGrid, yin=diag.cov, win=diag.cnt,xout=grid)
  }
  
  d <- length(grid)
  len <- grid[d] - grid[1]
  T.min <- min(which(grid >= grid[1] + 0.25 * len))
  T.max <- max(which(grid <= grid[d] - 0.25 * len))
  Diag <- (DiagCov - diag(CrCov))[T.min:T.max]
  w <- refund:::quadWeights(grid[T.min:T.max])
  sigma2 <- max((2/len) * drop(w %*% Diag), 1e-6)

  result <- list(CrCov=CrCov, sigma=sigma2)
  
  return(result)
}
