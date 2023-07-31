CompCrit = function (megaCov, a, W, M)
{
  J = length(megaCov)
  C = matrix(0, J, J)
  for (j in 1:J) {
    for (k in (1:J)[-j]) {
      projval <- megaCov[[j]][[k]] %*% W[[k]] %*% M[[k]] %*% a[[k]]
      C[j, k] <- drop(crossprod(M[[j]] %*% a[[j]], W[[j]] %*% projval))
    }
  }
  return(C)
}
