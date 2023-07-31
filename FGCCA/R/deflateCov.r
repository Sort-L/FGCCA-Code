deflateCov = function(megaCov, a, which.dfl, delfType, W)
{
  J <- length(megaCov)
  pjs <- sapply(a, length)

  deflatedMegaCov <- lapply(seq_len(J), function(j){list()})
  for (j in seq_len(J)){
    for (k in seq_len(J)) {
      deflCv <- megaCov[[j]][[k]]
      if (j %in% which.dfl) {
        if (delfType == "uncor"){
          normj <- t(a[[j]]) %*%  W[[j]] %*% megaCov[[j]][[j]] %*% W[[j]] %*% a[[j]]
          deflCv <- (diag(pjs[j]) - (megaCov[[j]][[j]] %*% W[[j]] %*% tcrossprod(a[[j]]) %*% W[[j]]) / as.numeric(normj)) %*% deflCv
        }
        if (delfType == "ortho") {
          deflCv <- (diag(pjs[j]) - (tcrossprod(a[[j]]) %*% W[[j]] / drop(crossprod(a[[j]], W[[j]] %*% a[[j]])))) %*% deflCv
        }
      }
      if (k %in% which.dfl) {
        if (delfType == "uncor"){
          normk <- t(a[[k]]) %*% W[[k]] %*% megaCov[[k]][[k]] %*% W[[k]] %*% a[[k]]
          deflCv <- deflCv %*% W[[k]] %*% (diag(pjs[k]) %*% solve(W[[k]]) - (tcrossprod(a[[k]]) %*% W[[k]] %*% megaCov[[k]][[k]]) / as.numeric(normk))
        }
        if (delfType == "ortho") {
          deflCv <- deflCv %*% W[[k]] %*% (diag(pjs[k]) %*% solve(W[[k]]) - (tcrossprod(a[[k]]) / drop(crossprod(a[[k]], W[[k]] %*% a[[k]]))))
        }
      }
      deflatedMegaCov[[j]][[k]] <- deflCv
    }
  }
  
  
  return(deflatedMegaCov)
}
