rm(list = ls())
PATH = paste0(Sys.getenv("PWD"), "/scripts/utils/")
library(foreach)

file.sources = list.files(path = PATH, pattern = "*.R", full.names = T)
sapply(file.sources, source, .GlobalEnv)

fgcca.method = function(data.obs, grids, M) {
  J <- length(grids)
  n <- dim(data.obs[[1]])[1]
  FPCAoptns <- lapply(seq_len(J), function(j){list(dataType="Sparse", nRegGrid=50)})
  
  Lys <- lapply(seq_len(J), function(j){lapply(seq_len(n), function(i){data.obs[[j]][i,]})})
  Lts <- lapply(seq_len(J), function(j){lapply(seq_len(n), function(i){grids[[j]]})})
  
  fgcca <- FGCCA::fgcca(Lys, Lts, ncomp=rep(M, J), verbose=F, scale=F, deflType="ortho",
                 optns=list(smoothing="lwls", scoreConnection=pracma::ones(J)))
  xi <- lapply(seq_len(M), function(m){sapply(seq_len(J), function(j){fgcca$Y[[j]][,m]})})
  phi <- lapply(seq_len(M), function(m){sapply(seq_len(J), function(j){fgcca$a[[j]][,m]})})
  data <- lapply(seq_len(J), function(j){t(fgcca$mus[[j]] + t(fgcca$Y[[j]] %*% t(fgcca$a[[j]])))})
  
  return(list(data=data))
}

mfpca.method = function(data.obs, grids, M) {
  J <- length(grids)
  n <- dim(data.obs[[1]])[1]
  
  funDatas <- lapply(seq_len(J), function(j){funData::funData(grids[[j]], data.obs[[j]])})
  mfunDatas <- funData::multiFunData(funDatas) 
  
  mfpca <- MFPCA::MFPCA(mfunDatas, M, uniExpansions = list(list(type="uFPCA"), list(type="uFPCA"), list(type="uFPCA")))
  
  data <- lapply(seq_len(J), function(j){t(drop(attr(mfpca$meanFunction[[j]], "X")) + t(mfpca$scores %*% attr(mfpca$functions[[j]], "X")))})

  return(list(data=data))
}

simulate = function(J, M, n, sigma, sparsity, methods) {
  grids <- lapply(seq_len(J), function(j){seq(0, 1, length.out=50)})
  
  sim <-  funData::simMultiFunData(type = "split", argvals = grids, M = M, eFunType = "Poly", eValType = "linear", N = n)
  simDatas <- lapply(seq_len(J), function(j){funData::sparsify(sim$simData[[j]], sparsity[1]*length(grids[[j]]), sparsity[2]*length(grids[[j]]))})
  data.true <- lapply(seq_len(J), function(j){t(sapply(seq_len(n), function(i){attr(sim$simData[[j]], "X")[i,]}))})
  data.true.na <- lapply(seq_len(J), function(j){t(sapply(seq_len(n), function(i){attr(simDatas[[j]], "X")[i,]}))})

  data.obs <- lapply(seq_len(J), function(j){data.true.na[[j]] + sapply(seq_len(length(grids[[j]])), function(i){rnorm(n, sd=sigma[j])})})
  apply.method = function(method, data.obs, grids, M) {
    tryCatch({
      results <- method(data.obs, grids, M)
      mrse <- MRSE(data.true, results$data, grids)
      return(mrse)
      }, 
      error=function(e) {
        cat("method did not work\n")
        return(NA)
      }
    )
    }
  
  mrse <- list()
  for (k in seq_len(length(methods))) {
	  metric <- apply.method(methods[[k]], data.obs, grids, M)
	  mrse[[names(methods)[k]]] <- metric
  }
  return(list(mrse = mrse))
}

methods <- list(FGCCA=fgcca.method, MFPCA=mfpca.method)
sparsity.vec <- c("Dense", "Low Sparsity", "Medium Sparsity", "High Sparsity")
n.vec <- c(100, 400, 900, 1600)

J <- 3
M <- 6
sigma <- 0.1
sigmas <- rep(sigma, J)
gridSizes <- rep(50, J)
N.SIM = 100

cat("Setting up cluster\n")
clst <- parallel::makeCluster(16)
doParallel::registerDoParallel(clst)
cat("Cluster setup\n")

results.mrse <- list()
for (sparsity in sparsity.vec) {
  results.mrse[[sparsity]] <- list()
  for (n in n.vec) {
    sparsity.v <- sparsityFromName(sparsity)

    runs <- foreach(i_sim=seq_len(N.SIM)) %dopar% {simulate(J, M, n, sigmas, sparsity.v, methods)}
    # runs <- lapply(seq_len(N.SIM), function(j){simulate(J, M, n, sigmas, sparsity.v, methods)})

    cat(paste0("Simulation : n = ", n, " ; sparsity = ", sparsity, "\n"))

    results.mrse[[sparsity]][[as.character(n)]] <- na.omit(data.frame(sapply(names(methods), function(method){
      sapply(runs, function(r){r$mrse[[method]]})}), n=n, M=M, sigma=sigma, sparsity=sparsity))
  }
}

parallel::stopCluster(clst)
cat("Cluster stopped")

tables.path <- paste0(PATH, "../../results/tables/design3/") 
figures.path <- paste0(PATH, "../../results/figures/design3/") 
filename <- paste0("sigma", format(sigma, nsmall=1))

save.results(results.mrse, NULL, NULL, tables.path, filename)
generate.figures(tables.path, figures.path, filename, "100", only.mrse=T)
