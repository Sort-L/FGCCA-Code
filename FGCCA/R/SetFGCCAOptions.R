SetFGCCAOptions <- function(Lys, Lts, optns){
  nbasis = optns$nbasis
  scoreMethod = optns$scoreMethod
  eigTreshold = optns$eigTreshold
  scoreConnection = optns$scoreConnection
  useReg = optns$useReg
  bws = optns$bws
  smoothing = optns$smoothing
  gridSize = optns$gridSize

  J <- length(Lys)
  which.fun <- which(!sapply(Lts, is.null))
  which.mtv <- which(sapply(Lts, is.null))

  if(is.null(nbasis)){
    nbasis = 10;
  }
  if(is.null(scoreMethod)){
    scoreMethod = rep("CE", length(Lys));
  } else {
    if (length(scoreMethod) == 1) scoreMethod = rep(scoreMethod, length(Lys));
  }
  if(is.null(eigTreshold)){
    eigTreshold = 0;
  }
  if(is.null(scoreConnection)){
    scoreConnection = ones(length(Lys));
  }
  if(is.null(useReg)){
    useReg = F;
  }
  if(is.null(bws)){
    bws = sapply(seq_len(J), function(j){ifelse(j %in% which.mtv, NA, 0.1*(max(unlist(Lts[[j]])) - min(unlist(Lts[[j]]))))})
  }
  if(is.null(smoothing)){
    smoothing = "gam"
  }
  if(is.null(gridSize)){
    gridSize = sapply(seq_len(J), function(j){ifelse(j %in% which.mtv, dim(Lys[[j]])[2], 50)})
  }
  
  retoptns <- list(nbasis = nbasis, scoreMethod = scoreMethod, eigTreshold = eigTreshold, scoreConnection = scoreConnection, 
                   useReg = useReg, bws=bws, smoothing=smoothing, gridSize=gridSize)
  
  invalidNames <- !names(optns) %in% names(retoptns)
  if (any(invalidNames)) {
    stop(sprintf('Invalid option names: %s',
                 paste0(names(optns)[invalidNames], collapse=', ')))
  }
  return(retoptns)
}
