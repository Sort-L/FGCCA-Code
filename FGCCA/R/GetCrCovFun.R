GetCrCovFun <- function(Ly1, 
                        Lt1, 
                        grid1, 
                        Ly2, 
                        Lt2, 
                        grid2, 
                        smoothing="gam",
                        nbasis=10,
                        bw1=NULL,
                        bw2=NULL) {
  obsGrid1 = sort(unique(unlist(Lt1)))
  obsGrid2 = sort(unique(unlist(Lt2)))
  
  N <- length(Ly1)
  p1 <- length(obsGrid1)
  p2 <- length(obsGrid2)
  
  cov.sum = cov.count = cov.mean = matrix(0, p1, p2)
  for (i in 1:N) {
    nas.points.1 = which(!is.na(Ly1[[i]]))
    nas.points.2 = which(!is.na(Ly2[[i]]))
    obs.points.1 = which(obsGrid1 %in% Lt1[[i]][nas.points.1])
    obs.points.2 = which(obsGrid2 %in% Lt2[[i]][nas.points.2])
    cov.count[obs.points.1, obs.points.2] = cov.count[obs.points.1, obs.points.2] + 1
    cov.sum[obs.points.1, obs.points.2] = cov.sum[obs.points.1, obs.points.2] +
      tcrossprod(Ly1[[i]][nas.points.1], Ly2[[i]][nas.points.2])
  }
  cov.mean = ifelse(cov.count == 0, NA, cov.sum/cov.count)
  
  row.vec = rep(obsGrid1, p2)
  col.vec = rep(obsGrid2, each=p1)
  cov.vec = as.vector(cov.mean)
  cnt.vec = as.vector(matrix(1, p1, p2))
  
  row.vec.pred = rep(grid1, length(grid2))
  col.vec.pred = rep(grid2, each=length(grid1))
  
  if (smoothing == "gam") {
    CrCov = matrix(predict(mgcv::gam(cov.vec ~ te(row.vec, col.vec, k=nbasis), weights = cnt.vec),
                           newdata = data.frame(row.vec = row.vec.pred, col.vec = col.vec.pred)), nrow=length(grid1), ncol=length(grid2))
  } else {
    nonzeros.vals = which(!is.na(cov.vec))
    pairs.df <- cbind(row.vec[nonzeros.vals], col.vec[nonzeros.vals])
    CrCov = fdapace::Lwls2D(c(bw1, bw2), kern="gauss", xin=pairs.df, yin=cov.vec[nonzeros.vals],
                            win = cnt.vec[nonzeros.vals], xout1=grid1, xout2=grid2, crosscov = T)
  }
  return(CrCov)
}
