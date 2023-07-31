rm(list = ls())
PATH = paste0(Sys.getenv("PWD"), "/scripts/utils/")

file.sources = list.files(path = PATH, pattern = "*.R", full.names = T)
sapply(file.sources, source, .GlobalEnv)
library(foreach)

fgcca.method.ce = function(data.obs, grids, M) {
  J <- length(grids)
  n <- dim(data.obs[[1]])[1]
  
  Lys <- lapply(seq_len(J), function(j){lapply(seq_len(n), function(i){data.obs[[j]][i,]})})
  Lts <- lapply(seq_len(J), function(j){lapply(seq_len(n), function(i){grids[[j]]})})
  
  fgcca <- FGCCA::fgcca(Lys, Lts, ncomp=rep(M, J), verbose=F, scale=F,
                 optns=list(smoothing="gam", scoreConnection=pracma::ones(J), scoreMethod="CE"))
  xi <- lapply(seq_len(M), function(m){sapply(seq_len(J), function(j){fgcca$Y[[j]][,m]})})
  phi <- lapply(seq_len(M), function(m){sapply(seq_len(J), function(j){fgcca$a[[j]][,m]})})
  data <- lapply(seq_len(J), function(j){t(fgcca$mus[[j]] + t(fgcca$Y[[j]] %*% t(fgcca$a[[j]])))})
  
  return(list(phi=phi, xi=xi, data=data))
}

fgcca.method.in = function(data.obs, grids, M) {
  J <- length(grids)
  n <- dim(data.obs[[1]])[1]
  
  Lys <- lapply(seq_len(J), function(j){lapply(seq_len(n), function(i){data.obs[[j]][i,]})})
  Lts <- lapply(seq_len(J), function(j){lapply(seq_len(n), function(i){grids[[j]]})})
  
  fgcca <- FGCCA::fgcca(Lys, Lts, ncomp=rep(M, J), verbose=F, scale=F,
                        optns=list(smoothing="gam", scoreConnection=pracma::ones(J), scoreMethod="IN"))
  xi <- lapply(seq_len(M), function(m){sapply(seq_len(J), function(j){fgcca$Y[[j]][,m]})})
  phi <- lapply(seq_len(M), function(m){sapply(seq_len(J), function(j){fgcca$a[[j]][,m]})})
  data <- lapply(seq_len(J), function(j){t(fgcca$mus[[j]] + t(fgcca$Y[[j]] %*% t(fgcca$a[[j]])))})
  
  return(list(phi=phi, xi=xi, data=data))
}

simulate = function(J, M, n, sigmas, sparsity, methods) {
  gridSizes <- rep(51, J)
  grids <- lapply(seq_len(J), function(j){seq(0, 1, length.out=gridSizes[j])})
  
  sim <-  simMultiblockFunData(argvals = grids, M = M, eFunType = "Fourier", eValType = "exponential", N = n)
  simDatas <- lapply(seq_len(J), function(j){funData::sparsify(sim$simData[[j]], sparsity[1]*length(grids[[j]]), sparsity[2]*length(grids[[j]]))})
  data.true <- lapply(seq_len(J), function(j){t(sapply(seq_len(n), function(i){attr(sim$simData[[j]], "X")[i,]}))})
  data.true.na <- lapply(seq_len(J), function(j){t(sapply(seq_len(n), function(i){attr(simDatas[[j]], "X")[i,]}))})
  
  xi.trues <- lapply(seq_len(M), function(m){sapply(seq_len(J), function(j){sim$trueScores[[j]][,m]})})
  phi.trues <- lapply(seq_len(M), function(m){sapply(seq_len(J), function(j){t(sim$trueFuns[[j]]@X)[,m]})})
  data.obs <- lapply(seq_len(J), function(j){data.true.na[[j]] + sapply(seq_len(length(grids[[j]])), function(i){rnorm(n, sd=sigmas[j])})})
  
  err.xi <- err.phi <- mrse <- list()
  for (k in seq_len(length(methods))) {
    method <- methods[[k]](data.obs, grids, M)
    err.phi[[names(methods)[k]]] <- sapply(seq_len(M), function(m){MSE.PHI(phi.trues[[m]], method$phi[[m]], grids)})
    err.xi[[names(methods)[k]]] <- sapply(seq_len(M), function(m){MSE.XI(xi.trues[[m]], method$xi[[m]])})
    mrse[[names(methods)[k]]] <- MRSE(data.true, method$data, grids)
  }
  
  return(list(err.xi=err.xi, err.phi=err.phi, mrse=mrse))
}

methods <- list(Integration=fgcca.method.in, Bayesian=fgcca.method.ce)
sparsity.vec <- c("Dense", "Low Sparsity", "Medium Sparsity", "High Sparsity")
n.vec <- c(100, 400)

J <- 2
M <- 6
sigma <- 1
sigmas <- rep(sigma, J)
N.SIM = 100

cat("Setting up cluster")
clst <- parallel::makeCluster(8)
doParallel::registerDoParallel(clst)
cat("Cluster set up")

results.mrse <- results.err.phi <- results.err.xi <- list()
for (sparsity in sparsity.vec) {
  results.mrse[[sparsity]] <- results.err.phi[[sparsity]] <- results.err.xi[[sparsity]] <- list()
  for (n in n.vec) {
    sparsity.v <- sparsityFromName(sparsity)

    runs <- foreach(i_sim=seq_len(N.SIM), .options.snow=opts) %dopar% {simulate(J, M, n, sigmas, sparsity.v, methods)}
    #runs <- pbmcapply::pbmclapply(seq_len(N.SIM), function(j){simulate(J, M, n, sigmas, sparsity.v, methods)}, mc.cores=2)
    #runs <- lapply(seq_len(N.SIM), function(j){simulate(J, M, n, sigmas, sparsity.v, methods)})
    
    cat(paste0("Simulation : n = ", n, " ; sparsity = ", sparsity, "\n"))
    
    results.mrse[[sparsity]][[as.character(n)]] <- data.frame(sapply(names(methods), function(method){
      sapply(runs, function(r){r$mrse[[method]]})}), n=n, M=M, sigma=sigma, sparsity=sparsity)
    results.err.phi[[sparsity]][[as.character(n)]] <- data.frame(sapply(names(methods), function(method){
      sapply(runs, function(r){r$err.phi[[method]]})}), n=n, M=M, sigma=sigma, sparsity=sparsity, m=rep(seq_len(M), N.SIM))
    results.err.xi[[sparsity]][[as.character(n)]] <- data.frame(sapply(names(methods), function(method){
      sapply(runs, function(r){r$err.xi[[method]]})}), n=n, M=M, sigma=sigma, sparsity=sparsity, m=rep(seq_len(M), N.SIM))
  }
}

parallel::stopCluster(clst)
cat("Cluster stopped")

tables.path <- paste0(PATH, "../../results/tables/design1/") 
figures.path <- paste0(PATH, "../../results/figures/design1/") 
filename <- paste0("sigma", format(sigma, nsmall=1))

save.results(results.mrse, results.err.phi, results.err.xi, tables.path, filename)
generate.figures(tables.path, figures.path, filename, "100")
