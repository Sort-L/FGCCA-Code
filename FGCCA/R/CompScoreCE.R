#' @export
CompScoreCE = function(Lys,
                       Lts,
                       a,
                       sigmas,
                       grids, 
                       megaCov,
                       W,
                       C_scores=ones(length(Lys)),
                       useReg=F,
                       scaleValues=rep(1, lengths(Lys))) {
  J <- length(Lys) # number of processes
  N <- length(Lys[[1]]) # number of individuals
  M <- NCOL(a[[1]]) # number of components
  pjs <- sapply(Lys, NCOL)
  
  which.fun <- which(sapply(seq_len(J), function(j){!is.matrix(Lys[[j]])}))
  which.mtv <- which(sapply(seq_len(J), function(j){is.matrix(Lys[[j]])}))
  lys <- Lys
  Lys <- lapply(seq_len(J), function(j){lapply(seq_len(N), function(i){if (is.matrix(Lys[[j]])) Lys[[j]][i,] else Lys[[j]][[i]]})})
  obs.points <- lapply(seq_len(J), function(j) lapply(seq_len(N), function(i){Lts[[j]][[i]][!is.na(Lys[[j]][[i]])]}))
  obs.values <- lapply(seq_len(J), function(j) lapply(seq_len(N), function(i){Lys[[j]][[i]][!is.na(Lys[[j]][[i]])]}))
  obsGrids <- lapply(seq_len(J), function(j){sort(unique(unlist(obs.points[[j]])))})
  
  for (j in which.fun) {
    if (max(obsGrids[[j]]) > max(grids[[j]])) warning(paste0("workGrid max (",  max(grids[[j]]), ") too small (", max(obsGrids[[j]]),"), removing values \n"))
    if (min(obsGrids[[j]]) < min(grids[[j]])) warning(paste0("workGrid min (",  min(grids[[j]]), ") too large (", min(obsGrids[[j]]),"), removing values \n"))
  }
  for (j in which.fun) {
    obs.values[[j]] <- lapply(seq_len(N), function(i){obs.values[[j]][[i]][obs.points[[j]][[i]] <= max(grids[[j]]) & obs.points[[j]][[i]] >= min(grids[[j]])]})
    obs.points[[j]] <- lapply(seq_len(N), function(i){obs.points[[j]][[i]][obs.points[[j]][[i]] <= max(grids[[j]]) & obs.points[[j]][[i]] >= min(grids[[j]])]})
  }
  
  obsGrids <- lapply(seq_len(J), function(j){sort(unique(unlist(obs.points[[j]])))})
  scores <- lapply(seq_len(J), function(j){matrix(0, N, M)})
  
  # Building mega covariance matrix on the observed time points and grid time points
  megaCovObs <- Sigmas <- lapply(seq_len(J), function(j){list()})
  StackSmoothCovLine <- muObs <- aObs <- StackSigmasLines <-  list()
  for (j in seq_len(J - 1)){
    for (k in seq(j + 1, J)){
      Sigmas[[j]][[k]] <- crossprod(a[[j]], W[[j]] %*% megaCov[[j]][[k]] %*% W[[k]] %*% a[[k]])
      Sigmas[[k]][[j]] <- t(Sigmas[[j]][[k]])
    }
  }
  
  for (j in which.fun) {
    aObs[[j]] <- fdapace::ConvertSupport(fromGrid = grids[[j]], toGrid=obsGrids[[j]], phi=a[[j]])
  }
  for (j in which.mtv) {
    aObs[[j]] <- a[[j]]
  }
  for (j in seq_len(J)) {
    Sigmas[[j]][[j]] <- crossprod(a[[j]], W[[j]] %*% (megaCov[[j]][[j]]) %*% W[[j]] %*% a[[j]])
    StackSigmasLines[[j]] <- do.call(cbind, Sigmas[[j]])
  }

  Sigma <- do.call(rbind, StackSigmasLines)
  eig <- eigen(Sigma)
  pos <- which(eig$values > 0)
  Sigma <- eig$vectors[,pos] %*% diag(eig$values[pos]) %*% t(eig$vectors[,pos])
  
  which.compute = which(!is.na(sigmas))
  
  # Computing the components
  # For each process
  for (j in which.fun) {
    blockIndices = which(C_scores[j,] == 1)
    # For each individual
    for (i in seq_len(N)){
      if (useReg | sigmas[j] == 0) {
        # if using regression approach
        procIdx <- which(obsGrids[[j]] %in% obs.points[[j]][[i]])
        sigma <- ifelse(is.na(sigmas[j]), 0, sigmas[j])
        Z.inv <- solve(crossprod(a[[j]]) + sigma * diag(M))
        scores[[j]][i, ] = Z.inv %*% t(aObs[[j]][procIdx,,drop=F]) %*% c(obs.values[[j]][[i]])
      } else {
        sigmaSelect <- unlist(lapply(blockIndices, function(k){(((k-1)*M+1):(k*M))}))
        SigmaBlock <- Sigma[(((j-1)*M+1):(j*M)),sigmaSelect]
        obsIndices <- lapply(seq_len(J), function(k){which(obsGrids[[k]] %in% obs.points[[k]][[i]])})
        Phi <- do.call(magic::adiag, lapply(blockIndices, function(k){aObs[[k]][obsIndices[[k]],,drop=F]}))
        sigmaDiag <- do.call(magic::adiag, lapply(blockIndices, function(k){sigmas[k]*diag(length(obs.points[[k]][[i]]))}))
        UiMui <- unlist(lapply(blockIndices, function(k){obs.values[[k]][[i]]}))
        scores[[j]][i, ] <- (1/scaleValues[j]) * SigmaBlock %*% t(Phi) %*% solve(Phi %*% Sigma[sigmaSelect, sigmaSelect] %*% t(Phi) + sigmaDiag) %*% UiMui
      }
    }
  }
  for (j in which.mtv) {
    scores[[j]] = lys[[j]] %*% a[[j]]
  }

  return(scores)
}
