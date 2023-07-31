MSE.XI = function(X_true, X_pred){
  sgn <- sign(drop(diag(crossprod(X_true, X_pred))))
  err <- mean(colMeans((X_true - t(sgn * t(X_pred)))**2))
  return(err)
}

MSE.PHI = function(X_true, X_pred, grids){
  sgn <- sign(drop(diag(crossprod(X_true, X_pred))))
  err <- (X_true - t(sgn * t(X_pred)))**2
  ierr <- mean(sapply(seq_len(length(grids)), function(j){fdapace::trapzRcpp(grids[[j]], err[,j])}))
  return(ierr)
}

MRSE = function(X_true, X_pred, grids){
  J <- length(X_pred)
  n <- dim(X_pred[[1]])[1]
  mrse <- (1/n) * sum(sapply(seq_len(n), function(i){sum(sapply(seq_len(J), function(j){fdapace::trapzRcpp(grids[[j]], (X_true[[j]][i,] - X_pred[[j]][i,])^2)}))}) / 
                        sapply(seq_len(n), function(i){sum(sapply(seq_len(J), function(j){fdapace::trapzRcpp(grids[[j]], (X_true[[j]][i,])^2)}))}))
  return(mrse)
}