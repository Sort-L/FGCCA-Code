CompScoreIN = function(Ly, Lt, a, grid, scaleValue) {
  N <- NROW(Ly) # number of individuals
  M <- NCOL(a)
  scores <- matrix(0, N, M)
  obsGrid <- sort(unique(unlist(Lt)))
  aObs <- fdapace::ConvertSupport(fromGrid=grid, toGrid=obsGrid, phi=a)
  for (i in seq_len(N)) {
    obs.points <- which(!is.na(Ly[[i]]))
    procIndices <- which(obsGrid %in% Lt[[i]][obs.points])
    scores[i,] <- apply(aObs[procIndices,,drop=F], 2, function(x){fdapace:::trapzRcpp(Lt[[i]][obs.points], 
                                                                               x * Ly[[i]][obs.points]) / scaleValue})
  }
  return(scores)
}
