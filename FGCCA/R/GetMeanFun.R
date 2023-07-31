GetMeanFun = function(Ly, Lt, 
                      grid,
                      smoothing="gam",
                      nbasis=10,
                      bw=NULL)
{
  obsGrid <- sort(unique(unlist(Lt)))
  N <- length(Ly)
  p <- length(obsGrid)
  
  val.vec <- unlist(Ly)
  val.sum <- val.count <- val.mean <- rep(0, p)
  
  for (i in seq_len(N)) {
    nas.points <- which(!is.na(Ly[[i]]))
    obs.points <- which(obsGrid %in% Lt[[i]][nas.points])
    val.count[obs.points] = val.count[obs.points] + 1
    val.sum[obs.points] = val.sum[obs.points] + Ly[[i]][nas.points]
  }
  val.mean <- ifelse(val.count == 0, NA, val.sum / val.count)
  
  obs.vec <- obsGrid
  val.vec <- as.vector(val.mean)
  cnt.vec <- as.vector(rep(1, p))
  
  if (smoothing == "gam") {
    mu <- as.vector(predict(mgcv::gam(val.vec ~ s(obs.vec, k=nbasis), weights = cnt.vec),
                            newdata=data.frame(obs.vec=grid)))
  }
  if (smoothing == "lwls") {
    nonzeros.vals = which(!is.na(val.vec))
    mu <- fdapace::Lwls1D(bw=bw, kern="gauss", win=cnt.vec[nonzeros.vals], 
                          xin=obs.vec[nonzeros.vals], yin=val.vec[nonzeros.vals], xout=grid)
  }
  if (smoothing == "ns") {
    mu <- colMeans(matrix(unlist(Ly), nrow=N, byrow=T), na.rm=T)
  }
  
  return(mu)
}