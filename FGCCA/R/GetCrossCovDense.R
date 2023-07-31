# <!> Function directly inspired from fdapace package
# This function obtains the sample covariance matrix at observed grid
# for dense regular functional data

######
# Input:
######  
#  ymat: n by p matrix of dense regular functional data
#  mu: p-dim vector, estimated cross-sectional mean
#  optns: options for FPCA function
#  y: list of amplitude information
#  t: list of time information
######
# Output: 
######
#  a SmoothCov object containing: 
#    - p by p matrix of sample cov surface estimation on observed grid
#    - NULL for all other entires
##########################################################################

GetCrossCovDense <- function(ymat1, mu1, ymat2, mu2, optns1, optns2) {
  if(!(optns1$dataType %in% c('Dense', 'DenseWithMV') & optns2$dataType %in% c('Dense', 'DenseWithMV'))){
    stop('Sample Covariance is only applicable for option: dataType = "Dense" or "DenseWithMV"!')
  }
  # if( optns$muCovEstMethod == 'cross-sectional' ){
  n = nrow(ymat1)
  m1 = ncol(ymat1)
  m2 = ncol(ymat2)
  if( !is.null(optns1$userMu) & !is.null(optns2$userMu) ) {
    ymat1 = ymat1 - matrix(rep(times= nrow(ymat1), mu1), ncol= ncol(ymat1), byrow=TRUE)
    ymat2 = ymat2 - matrix(rep(times= nrow(ymat2), mu2), ncol= ncol(ymat2), byrow=TRUE)
    K = matrix( rep(0,m1*m2), m1)
    for( i in (1:m1)){
      for( j in (1:m2)){
        XcNaNindx = which(is.na(ymat1[,i]));
        YcNaNindx = which(is.na(ymat2[,j]));
        NaNrows = union(XcNaNindx, YcNaNindx);
        # Find inner product of the columns with no NaN values
        indx = setdiff( 1:n, NaNrows)
        K[i,j] =  sum(ymat1[indx,i] * ymat2[indx,j]) * (1/(n-1-length(NaNrows)));  
      }
    }    
  } else {
    K = cov(ymat1, ymat2, use = 'pairwise.complete.obs') # sample variance using non-missing data
  }
  
  if (any(is.na(K))) {
    stop("Data is too sparse to be considered DenseWithMV. Remove sparse observations or specify dataType='Sparse' for FPCA") 
  }
  
  ret = list('rawCov' = NULL, 'smoothCrossCov' = K, 'bwCov' = NULL, outGrid = NULL)
  class(ret) = "SmoothCrossCov"
  return(ret)
  # } else {
  # stop('optns$muCovEstMethod is unknown!\n')
  # }
}