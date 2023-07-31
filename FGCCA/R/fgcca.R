#' Functional Generalized Canonical Correlation Analysis (FGCCA)
#' 
#' Functional Generalized Canonical Correlation Analysis (FGCCA) is a flexible framework for exploring relationships between any number J of random processes.
#' It is robust to sparse and irregular functional data.
#'
#' @param Lys A J list of n lists containing the observations
#' @param Lts A J list of n lists containing the time-points
#' @param connection The design matrix of the FGCCA model (J*J)
#' @param tau The regularization parameter for each process
#' @param scheme The scheme function (should be 'horst', 'factorial' or 'centroid')
#' @param init The initialization method (should be 'random' or 'svd')
#' @param tol Tolerance parameter for the algorithm
#' @param verbose A logical value indicating if showing computation details
#' @param center A logical value indicating if applying centering on the processes
#' @param scale A logical value indicating if applying a variance scaling on the processes
#' @param deflType The deflation type (should be 'uncor' (default) or 'ortho')
#' @param ncomp The number of components to extract for each process
#' @param optns A list of options
#'
#' @details Options for FGCCA are :
#' \describe{
#' Smoothing and score computation options (optns)
#' \item{smoothing}{"gam", "lwls" or "ns" (no smoothing), must be on a regular grid. Default is "gam".}
#' \item{gridSize}{Size of the smoothing grids. Default is 50.}
#' \item{scoreConnection}{Connection matrix used in the Gaussian conditioning. Default is full of 1 matrix. Can be useful to change for supervised learning}
#' \item{scoreMethod}{"CE" or "IN". "CE" is Conditional Expectation (Bayesian Approach). "IN" is by numerical INtegration Default is "CE".}
#' }
#'
#' @return A fgcca object containing the canonical functions (a) and the canonical components (Y)
#' @export
fgcca = function(Lys, Lts,
                 connection = 1 - diag(length(Lys)),
                 tau = rep(1, length(Lys)),
                 scheme = "factorial",
                 init = "random", 
                 tol = 1e-14,
                 deflType = "ortho", 
                 verbose = TRUE,
                 center = TRUE,
                 scale = FALSE,
                 ncomp = rep(1, length(Lys)),
                 optns = list())
{
  J <- length(Lys)
  N <- NROW(Lys[[1]])
  pjs <- sapply(Lys, NCOL)
  lys <- lts <- list()
  
  # Checks
  optns <- SetFGCCAOptions(Lys, Lts, optns=optns)
  match.arg(init, c("svd", "random"))
  match.arg(deflType, c("ortho", "uncor"))
  match.arg(optns$smoothing, c("gam", "lwls", "ns"))
  
  which.fun <- which(!sapply(Lts, is.null))
  which.mtv <- which(sapply(Lts, is.null))
  
  ndefl <- ncomp - 1
  ndeflMax <- max(ndefl)
  
  obsGrids <- workGrids <- W <- list()
  for (j in which.fun) {
    obsGrids[[j]] <- sort(unique(unlist(Lts[[j]])))
    if (optns$smoothing == "ns") {
      workGrids[[j]] <- obsGrids[[j]]
    } else {
      workGrids[[j]] <- seq(min(obsGrids[[j]]), max(obsGrids[[j]]), length.out=optns$gridSize[j])
    }
    W[[j]] <- diag(refund:::quadWeights(workGrids[[j]]))
    lts[[j]] <- Lts[[j]]
  }
  
  mus <- list()
  if (center) {
    for (j in which.fun) {
      mus[[j]] = GetMeanFun(Lys[[j]], Lts[[j]], workGrids[[j]], optns$smoothing, optns$nbasis, optns$bws[j])
      lys[[j]] <- lapply(seq_len(N), function(i){Lys[[j]][[i]] - pracma::interp1(workGrids[[j]], mus[[j]], Lts[[j]][[i]])})
    }
  }
  
  # Computing mean for multivariate blocks
  for (j in which.mtv) {
    obsGrids[[j]] <- as.double(1:pjs[j])
    workGrids[[j]] <- as.double(1:pjs[j])
    W[[j]] <- diag(length(workGrids[[j]]))
    mus[[j]] = colMeans(Lys[[j]])
    lys[[j]] <- scale(Lys[[j]], scale=T)
    lts[[j]] <- lapply(seq_len(N), function(i){as.double(1:pjs[j])})
  }
  
  if (verbose) {
    cat("Mus and grids computed\n")
  }
  
  if (verbose) {
    cat("Computation of the covariance and cross-covariance matrices \n")
  }
  
  megaCov <- lapply(1:J, function(k){list()})
  sigmas <- rep(0, J)
  for (j in seq_len(J)) {
    if (j %in% which.fun) {
      if (optns$smoothing != "ns") {
        C <- GetCovFun(lys[[j]], Lts[[j]], workGrids[[j]], optns$smoothing, optns$nbasis, optns$bws[j])
        megaCov[[j]][[j]] <- C$CrCov
        sigmas[j] <- C$sigma
      } else {
        ymat = fdapace:::List2Mat(lys[[j]], lts[[j]])
        Cov = cov(ymat, use = 'pairwise.complete.obs')
        ord <- 2 
        sigma2 <- mean(diff(t(ymat), differences=ord)^2, na.rm=TRUE) / choose(2 * ord, ord)
        megaCov[[j]][[j]] <- Cov
        sigmas[j] <- sigma2
      }
    }
    if (j %in% which.mtv) {
      megaCov[[j]][[j]] <- var(lys[[j]])
    }
  }
  
  scaleValues <- rep(1, J)
  if (scale) {
    for (j in which.fun) {
      scaleValues[j] <- sqrt(1/trapz(workGrids[[j]], diag(megaCov[[j]][[j]])))
      lys[[j]] <- lapply(lys[[j]], function(ly){scaleValues[j]*ly})
      megaCov[[j]][[j]] <- scaleValues[j] * scaleValues[j] * megaCov[[j]][[j]]
      sigmas[j] <- scaleValues[j] * scaleValues[j] * sigmas[j]
    }
  }
  
  for (j in seq_len(J - 1)) {
    for (k in seq(j + 1, J)) {
      if (j %in% which.fun && k %in% which.fun) {
        if (optns$smoothing != "ns") {
          CrCov <- GetCrCovFun(lys[[j]], lts[[j]], workGrids[[j]], lys[[k]], lts[[k]], workGrids[[k]], 
                               optns$smoothing, optns$nbasis, optns$bws[j], optns$bws[k])
        } else {
          y1mat = fdapace:::List2Mat(lys[[j]], lts[[j]])
          y2mat = fdapace:::List2Mat(lys[[k]], lts[[k]])
          CrCov = cov(y1mat, y2mat, use = 'pairwise.complete.obs')
        }
      }
      if (j %in% which.mtv && k %in% which.mtv) {
        CrCov <- cov(lys[[j]], lys[[k]])
      } 
      if (j %in% which.mtv && k %in% which.fun) {
        CrCov = t(GetMixCovMatrix(Ly.fun=lys[[k]], Lt.fun=lts[[k]], grid=workGrids[[k]], Ly.mtv=lys[[j]],
                                  optns$smoothing, optns$nbasis, optns$bws[k]))
      }
      if (j %in% which.fun && k %in% which.mtv) {
        CrCov = t(GetMixCovMatrix(Ly.fun=lys[[j]], Lt.fun=lts[[j]], grid=workGrids[[j]], Ly.mtv=lys[[k]],
                                  optns$smoothing, optns$nbasis, optns$bws[j]))
      }
      megaCov[[j]][[k]] <- CrCov #/ sqrt(length(workGrids[[j]])*length(workGrids[[k]]))
      megaCov[[k]][[j]] <- t(CrCov) #/ sqrt(length(workGrids[[j]])*length(workGrids[[k]]))
      }
    }
    
  pjs <- sapply(workGrids, length)
  
  ##### FGCCA Core : /!\ Based on "old" version of RGCCA 2.1.2 #####
  
  if (mode(scheme) != "function") {
    if (verbose)
      cat("Computation of the CovGCCA block components based on the",
          scheme, "scheme \n")
  }
  if (mode(scheme) == "function" & verbose)
    cat("Computation of the CovGCCA block components based on the g scheme \n")

  P <- a <- Y <- astar <- NULL
  deflMegaCov <- crit <- list()
  deflMegaCov[[1]] <- megaCov
  
  for (b in seq_len(J)) a[[b]]                     <- matrix(NA, pjs[b], ndeflMax + 1)
  for (b in seq_len(J)) astar[[b]]                 <- matrix(NA, pjs[b], ndeflMax + 1)
  for (b in seq_len(J)) Y[[b]]                     <- matrix(NA, N, ndeflMax + 1)
  for (b in seq_len(J)) deflMegaCov[[1]][[b]][[b]] <- deflMegaCov[[1]][[b]][[b]] + sigmas[b] * diag(pjs[b])

  # First order function estimation
  fgcca.result <- covgccak(megaCov, workGrids, connection, W, tau = tau, 
                           scheme = scheme, init = init, tol = tol, verbose = verbose)

  for (b in seq_len(J)) a[[b]][, 1]     <- fgcca.result$a[[b]]
  for (b in seq_len(J)) astar[[b]][, 1] <- fgcca.result$a[[b]]
  crit[[1]]                             <- fgcca.result$crit

  # Higher order function estimation
  if (ndeflMax > 0) {
    R <- Lys
    for (n in seq(2, (ndeflMax + 1))) {
      if (verbose)
        cat(paste0("Computation of the CovGCCA block components #", n, " is under
                 progress...\n"))
      which.dfl <- which(sapply(ncomp, function(nc){nc >= n}))
      deflMegaCov[[n]] <- deflateCov(deflMegaCov[[n-1]], lapply(a, function(x){x[, n-1, drop=FALSE]}), which.dfl, deflType, W)
      
      fgcca.result <- covgccak(deflMegaCov[[n]], workGrids, connection, W, tau = tau, 
                               scheme = scheme, init = init, tol = tol, verbose = verbose)

      crit[[n]] <- fgcca.result$crit

      for (b in 1:J) a[[b]][, n] <- fgcca.result$a[[b]]
      for (b in 1:J) 
      {
        if (deflType == "uncor") {
          astar[[b]][, n] <- (diag(pjs[b]) - tcrossprod(a[[b]][,n-1]) %*% W[[b]] %*% deflMegaCov[[n-1]][[b]][[b]] %*% W[[b]] / 
                                drop(crossprod(a[[b]][,n-1], W[[b]] %*% deflMegaCov[[n-1]][[b]][[b]] %*% W[[b]] %*% a[[b]][,n-1]))) %*% a[[b]][, n]
        } else {
          astar[[b]][, n] <- a[[b]][, n]
        }
      }
    }
  }
  
  if (verbose) 
    cat("Estimating components")
  # Component estimation
  if (any(optns$scoreMethod == "CE")) GPACE.results <- CompScoreCE(lys, lts, a, sigmas, workGrids, megaCov, W, optns$scoreConnection, optns$useReg, scaleValues)
  for (j in seq_len(J)) {
    if (optns$scoreMethod[j] == "CE") Y[[j]] <- GPACE.results[[j]]
    if (optns$scoreMethod[j] == "IN") Y[[j]] <- if (j %in% which.fun) CompScoreIN(lys[[j]], Lts[[j]], astar[[j]], workGrids[[j]], scaleValues[j]) else Lys[[j]] %*% a[[j]]
  }
  
  if (verbose) 
    cat("Deflating components")
  # Decorrelating components (for uncorrelated components (type 1) deflation type)
  Ystar <- Y
  if (deflType == "uncor")
    for (j in 1:J)
      for (n in seq(2, (ndeflMax + 1)))
        for (np in  seq(1, n-1))
          Ystar[[j]][,n] <- Ystar[[j]][,n] - (tcrossprod(Ystar[[j]][,np]) / drop(crossprod(Ystar[[j]][,np]))) %*% Y[[j]][,n]
  
  Y <- Ystar

  if (N == 0) {
    crit = unlist(crit)
  }

  out <- list(a = a, Y = Y, mus = mus, astar = astar, sigmas=sigmas, megaCov = megaCov, W = W, scaleValues = scaleValues,
              deflMegaCov = deflMegaCov, grids = workGrids, crit = crit, optns = optns)

  class(out) <- "fgcca"

  return(out)
}
