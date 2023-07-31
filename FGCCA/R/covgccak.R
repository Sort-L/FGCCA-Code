#' 
#' 
#' The function covgccak() is called by fgcca() and does not have to be used by
#' the user. The function returns the canonical functions/vectors for random object
#'
#' @param megaCov A list of lists containing the covariance/cross-covariance matrices
#' @param grids The work grid for each process
#' @param C The design matrix of the FGCCA model (J*J)
#' @param W A list of weighting matrices, allowing for functional data to perfom numerical integration
#' @param tau The regularization parameter for each process
#' @param scheme The scheme function (should be 'horst', 'factorial' or 'centroid')
#' @param verbose A logical value indicating if verbose
#' @param init The initialization method (should be 'random' or 'svd')
#' @param tol Tolerance parameter for the algorithm
#'
#' @return fgcca object#'
covgccak=function(megaCov, grids, C, W,
                  tau = rep(1, length(megaCov)),
                  scheme = "centroid", 
                  verbose = FALSE,
                  init = "random", 
                  tol = 1e-16)
{
  
  if(mode(scheme) != "function")
  {
    if(!scheme %in% c("horst", "factorial", "centroid"))
    {stop_rgcca("Please choose scheme as 'horst', 'factorial', 'centroid'")}
    if(scheme == "horst"){ g <- function(x) x }
    if(scheme == "factorial"){ g <- function(x)  x^2 }
    if(scheme == "centroid"){ g <- function(x) abs(x) }
  }
  else g <- scheme
  
  J <- length(megaCov) # number of blocks
  pjs <- sapply(grids, length) # number of time points per grid
  
  a <- Y <- M <- list()
  
  for (j in seq_len(J)) {
    if (init == "random") a[[j]] <- rnorm(pjs[j]) # random initialization
    if (init == "svd") a[[j]] <- Re(eigen(megaCov[[j]][[j]])$vectors[, 1]) # svd initialization
  }
  
  for (j in seq_len(J)) {
    a[[j]] <- (1 / sqrt(drop(crossprod(a[[j]], W[[j]] %*% a[[j]])))) * a[[j]]
    if (tau[j] == 1) {
      M[[j]] <- diag(pjs[j])
    } else {
      M[[j]] <- sqrtm(tau[j] * diag(pjs[j]) + ((1 - tau[j])) * megaCov[[j]][[j]])$Binv
      # a[[j]] <- drop(1/sqrt(t(a[[j]]) %*% t(W[[j]]) %*% M[[j]] %*% W[[j]] %*% a[[j]]))*(M[[j]] %*% W[[j]] %*% a[[j]])
    }
  }
  
  crit_old <- sum(C * g(CompCrit(megaCov, a, W, M)))
  iter = 1
  crit = numeric()
  a_old = a
  
  dg = Deriv::Deriv(g, env = parent.frame())
  
  repeat
  {
    for (j in 1:J) {
      grad <- matrix(0, pjs[j], 1)
      for (k in (1:J)[-j]) {
        proj <- megaCov[[j]][[k]] %*% W[[k]] %*% M[[k]] %*% a[[k]]
        cov <- dg(drop(crossprod(M[[j]] %*% a[[j]], W[[j]] %*% proj)))
        grad <- grad  + 2 * as.numeric(C[j, k] * cov) * proj
      }
      # if (tau[j] == 1) {
      a[[j]] <- grad / sqrt(drop(crossprod(grad, W[[j]] %*% grad)))
      # } else {
      #   a[[j]] <- drop(1 / sqrt(drop(crossprod(grad, t(W[[j]]) %*% M[[j]] %*% W[[j]] %*% grad)))) * (M[[j]] %*% W[[j]] %*% grad)
      # }
    }
    
    crit[iter] <- sum(C * g(CompCrit(megaCov, a, W, M)))
    
    if (verbose & (iter%%1) == 0)
    {
      cat(" Iter: ", formatC(iter, width = 3, format = "d"),
          " Fit:", formatC(crit[iter], digits = 8,
                           width = 10, format = "f"),
          " Dif: ", formatC(crit[iter] - crit_old, digits = 8,
                            width = 10, format = "f"), "\n")
    }
    stopping_criteria = c(drop(crossprod(Reduce("c", mapply("-", a, a_old)))),
                          abs(crit[iter] - crit_old))
    
    if (any(stopping_criteria < tol) | (iter > 1000))
    {break}
    crit_old = crit[iter]
    a_old <- a
    iter <- iter + 1
  }
  
  if (iter > 1000)
    warning("The CovGCCA algorithm did not converge after 1000 iterations.")
  if (iter < 1000 & verbose)
    cat("The CovGCCA algorithm converged to a stationary point after",
        iter - 1, "iterations \n")
  if (verbose)
  {
    plot(crit[1:iter], xlab = "iteration", ylab = "criteria")
  }
  
  a <- lapply(1:J, function(j){M[[j]] %*% a[[j]]})
  
  result <- list(a = a, Y = Y, crit = crit)
  return(result)
}
