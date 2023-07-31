simMultiblockFunData <- function(argvals, M, eFunType, ignoreDeg = NULL, eValType, N)
{
  if(! (is.list(argvals) & all(is.numeric(unlist(argvals)))) )
    stop("Parameter 'argvals' must be passed as a list of numerics.")
  
  if(! all(is.numeric(unlist(M))))
    stop("Parameter 'M' must contain only numerics.") 
  
  if(! all(is.character(unlist(eFunType))))
    stop("Parameter 'eFunType' must contain only strings.")
  
  if(!(is.null(ignoreDeg ) | all(is.numeric(ignoreDeg), ignoreDeg > 0)))
    stop("Parameter 'ignoreDeg' must be either NULL or a vector of positive numbers.") 
  
  if(! all(is.character(eValType), length(eValType) == 1))
    stop("Parameter 'eValType' must be passed as a string.")
  
  if(! all(is.numeric(N), length(N) == 1, N > 0))
    stop("Parameter 'N' must be passed as a positive number.") 
  
  # number of eigenfunctions generated
  Mtotal <- M
  
  # number of elements in multivariate functional basis
  p <- length(argvals)
  
  # generate eigenfunctions
  trueFuns <- lapply(seq_len(p), function(j){funData::eFun(argvals[[j]], M, ignoreDeg, eFunType)})
  
  # generate eigenvalues and scores
  trueVals <- funData:::eVal(Mtotal, eValType)
  cov <- diag(rep(trueVals, p))
  for (j in seq_len(p)) for (k in seq_len(p)) if (j != k) cov[((j-1)*Mtotal+1):(j*Mtotal), ((k-1)*Mtotal+1):(k*Mtotal)] = diag(trueVals) / 1.5
  scoresjoint <- MASS::mvrnorm(N, mu=rep(0, p*Mtotal), Sigma = cov)
  scores <- lapply(seq_len(p), function(j){scoresjoint[,((j-1)*Mtotal+1):(j*Mtotal)]})
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    X <- apply(trueFuns[[j]]@X, -1, function(v){scores[[j]] %*% v})
    
    if(N == 1)
      dim(X) <- c(1, funData:::nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData:::funData(trueFuns[[j]]@argvals, X)
  } 
  
  return(list(simData = funData:::multiFunData(simData),
              trueFuns = trueFuns,
              trueVals = trueVals,
              trueScores = scores))
}

sparsityFromName = function(name) {
  if (name == "High Sparsity")  return(c(0.1, 0.4))
  if (name == "Medium Sparsity")  return(c(0.4, 0.8))
  if (name == "Low Sparsity")  return(c(0.8, 1.0))
  if (name == "Dense")    return(c(1.0, 1.0))
}