#' Forecasting h-steps ahead using Recursive Least Squares Fast
#'
#' Consider the following LS-fitted Model with intercept:
#' y_(t+h) = beta_0 + x_(jt) * beta + u_(t+h)
#' which is used to generate out-of-sample forecasts of y, h-steps ahead (h=1,2,3,. . . ).
#' Notes: (1) first estimation window is (1,...,k0) and last window is 
#' (1,....,n-h) for k0 = round(n*pi0). First forecast is yhat(k0+h|k0)
#' and last forecast is yhat(n|n-h). There are a total of (n-h-k0+1)
#' forecasts and corresponding forecast errors. (2) this fast version of the
#' recursive least squares algorithm uses the Sherman-Morrison matrix
#' formula to avoid matrix inversions at each recursion. (3) x_(jt) is the j^th predictor in x (j^th column).
#'
#' recursive_hstep_fast is the fast version that avoids the recursive calculation of inverse of the matrix using Sherman-Morrison formula.
#'
#' @param y an outcome series, which should be numeric and one dimensional.
#' @param x a predictor matrix (intercept would be added automatically).
#' @param pi0 Fraction of the sample, which should be within 0 and 1.
#' @param h Number of steps ahead to predict, which should be a positive integer.
#' @return Series of residuals estimated
#' @examples
#' x<- rnorm(15);
#' y<- x+rnorm(15);
#' temp1 <- recursive_hstep_fast(y,x,pi0=0.5,h=1);
#' @export



recursive_hstep_fast=function(y,x,pi0,h)
{
  #Stopping criteria
  if(is.null(y)){
    stop("y must be one dimension")
  }
  
  y <- as.matrix(y)
  ny=dim(y)[1]
  py=dim(y)[2]
  
  if(py > 1){
    stop("y must be one dimension")
  }
  
  if(anyNA(as.numeric(y))){
    stop("y must not contain NA")
  }
  
  if(is.null(x)){
    stop("x must be one dimension")
  }
  
  x <- as.matrix(x)
  n=dim(x)[1]
  p=dim(x)[2] #[n,p] = size(X);
  
  if(anyNA(as.numeric(x))){
    stop("x must not contain NA")
  }
  
  if(ny != n){
    stop("y and x must have same length of rows")
  }
  
  if(is.null(pi0)){
    stop("pi0 must be between 0 and 1")
  }else if(!is.numeric(pi0)){
    stop("pi0 must be between 0 and 1")
  }
  
  if(pi0 >= 1 || pi0 <= 0){
    stop("pi0 must be between 0 and 1")
  }
  
  if(is.null(h)){
    stop("h must be a positive integer")
  }else if(!is.numeric(h)){
    stop("h must be a positive integer")
  }
  
  h <- as.integer(h)
  
  if(h <= 0 || h > (n-1)){
    stop("h must be a positive integer")
  }
  
  
  #Initialisation of parameters
  xx=as.matrix(x)
  k0 = round(n * pi0)
  ehat = matrix(NA,n-k0-h+1,p+1) #
  
  #for ehat0 and ehatj
  for (j in 1:(p+1)){
    
    if (j < (p+1)){
      x= cbind(rep(1,n),xx[,j])
    }else{
      x=cbind(rep(1,n))
    }
    
    nc=dim(x)[2]
    M = array(NA, dim = c(nc,nc,n-k0+1))
    beta = matrix(NA,nc,n-k0+1) #
    
    
    #Calculate the first M and Beta at k0
    iXmatk0 = solve(t(x[1:(k0-h),])%*%x[1:(k0-h),]) #M_k0
    bhat_k0 = iXmatk0%*%(t(x[1:(k0-h),])%*%y[(1+h):(k0)]) #Beta_k0
    
    M[,,1] = iXmatk0
    beta[,1] = bhat_k0
    
    
    #iteractively update M and Beta for (k0+h) to n
    for (t in 1:(n-k0)){
      M[,,t+1] = M[,,t]-((M[,,t]%*%x[t+k0-h,]%*%t(x[t+k0-h,])%*%M[,,t])/c(1+t(x[t+k0-h,])%*%M[,,t]%*%x[t+k0-h,]))
      beta[,t+1] = beta[,t]+(M[,,t]%*%x[t+k0-h,]%*%(y[t+k0]-x[t+k0-h,]%*%beta[,t]))/c(1+t(x[t+k0-h,])%*%M[,,t]%*%x[t+k0-h,])
    }
    
    #update ehat from k0+h to n
    for (s in (k0+h):n){
      ehat[s-k0-h+1,j] = y[s]-x[s-h,]%*%beta[,s-k0-h+1]
    }
    
  }
  
  ehat0 = ehat[,p+1]
  ehatj = ehat[,1:p]

  
  return(list(ehat0=ehat0,ehatj=ehatj))
}
