#' Constructs the DMBAR Test statistic in GP2023
#'
#' Consider a linear predictive regression setting with a potentially large set of candidate predictors.
#' This work is concerned with detecting the presence of out of sample predictability based on out of sample MSE comparisons.
#' For details of the test, please refer to Gonzalo and Pitarakis (2023).
#'
#' @param ehat0 n by 1 vector of out of sample forecast errors from benchmark model with only intercept.
#' @param ehatj n by j vector of out of sample forecast errors from models (j=1,...,p) estimated with one predictor (j) per time.
#' @param mu0 sample split parameter (must be different from 0.5).
#' @param pvcutoffkp pvalue cutoff used to decide whether the global null is rejected when identifying the key player conditional on rejecting the global null.
#' @return A list of Test statistic, pvalue and key player across 4 alternative formulations of the test statistics (lrvar under 0 vs 1; power enhanced vs non-power enhanced (notation: 0, 1, 0adj, 1adj).
#' @references Gonzalo, J., & Pitarakis, J. Y. (2023). Out-of-sample predictability in predictive regressions with many predictor candidates. International Journal of Forecasting, 1166-1178.
#' @import stats
#' @examples
#' ehat0<- rnorm(15);
#' ehatj<- rnorm(15);
#' temp <- DMBAR_Test(ehat0,ehatj,mu0=0.4,pvcutoffkp=0.1);
#' @export


DMBAR_Test=function(ehat0, ehatj, mu0, pvcutoffkp)
{

  #Stopping criteria
  if(is.null(ehat0)){
    stop("ehat0 must be one dimension")
  }

  ehat0 <- as.matrix(ehat0)
  nehat0=dim(ehat0)[1]
  pehat0=dim(ehat0)[2]

  if(pehat0 > 1){
    stop("ehat0 must be one dimension")
  }

  if(anyNA(as.numeric(ehat0))){
    stop("y must not contain NA")
  }

  if(is.null(ehatj)){
    stop("ehatj must be at least one dimension")
  }

  ehatj <- as.matrix(ehatj)
  nehatj=dim(ehatj)[1]
  pehatj=dim(ehatj)[2] #[n,p] = size(X);

  if(anyNA(as.numeric(ehatj))){
    stop("ehatj must not contain NA")
  }

  if(nehat0 != nehatj){
    stop("ehat0 and ehatj must have same length of rows")
  }

  if(is.null(mu0)){
    stop("mu0 must be between 0 and 1 but not 0.5")
  }else if(!is.numeric(mu0)){
    stop("mu0 must be between 0 and 1 but not 0.5")
  }

  if(mu0 >= 1 || mu0 <= 0 || mu0 == 0.5){
    stop("mu0 must be between 0 and 1 but not 0.5")
  }

  if(is.null(pvcutoffkp)){
    stop("pvcutoffkp must be between 0 and 1")
  }else if(!is.numeric(pvcutoffkp)){
    stop("pvcutoffkp must be between 0 and 1")
  }

  if(pvcutoffkp >= 1 || pvcutoffkp <= 0 ){
    stop("pvcutoffkp must be between 0 and 1")
  }

  #Initialisation of parameters

  ns = nehatj
  p = pehatj
  m0 = round(mu0*ns)
  kappa0 = ((1-2*mu0)^2)/(4*mu0*(1-mu0))

  ehat0sq = ehat0^2
  ehatjsq = ehatj^2
  ehatjsqadj = ehatjsq - (ehat0-ehatj)^2

  phihatsq0 = t(ehat0sq - mean(ehat0sq))%*%(ehat0sq - mean(ehat0sq))/ns
  omhatsq0 = phihatsq0 * kappa0  # lrvar using null residuals

  #Preallocations
  phihatsqj = matrix(NA,p,1)
  omhatsqj = matrix(NA,p,1)

  gpstat0 = matrix(NA,p,1)
  gpstat1 = matrix(NA,p,1)

  gpstat0adj = matrix(NA,p,1)
  gpstat1adj = matrix(NA,p,1)


  for (s in 1:p) {
  phihatsqj[s,] = t(ehatjsq[,s] - mean(ehatjsq[,s]))%*%(ehatjsq[,s] - mean(ehatjsq[,s]))/ns
  omhatsqj[s,] = phihatsqj[s,]*kappa0 # lrvars using residuals under alternative
  }


  #test statistics (differ depending on lrvar choice 0/1 and whether adjusted or not)

  for (s in 1:p){
  gpstat0[s] = (sqrt(ns)*(0.5*(mean(ehat0sq[1:m0]) + mean(ehat0sq[(m0+1):ns])) - mean(ehatjsq[,s])))/sqrt(omhatsq0)
  gpstat1[s] = (sqrt(ns)*(0.5*(mean(ehat0sq[1:m0]) + mean(ehat0sq[(m0+1):ns])) - mean(ehatjsq[,s])))/sqrt(omhatsqj[s])
  gpstat0adj[s] = (sqrt(ns)*(0.5*(mean(ehat0sq[1:m0]) + mean(ehat0sq[(m0+1):ns])) - mean(ehatjsqadj[,s])))/sqrt(omhatsq0)
  gpstat1adj[s] = (sqrt(ns)*(0.5*(mean(ehat0sq[1:m0]) + mean(ehat0sq[(m0+1):ns])) - mean(ehatjsqadj[,s])))/sqrt(omhatsqj[s])
  }

  #Aggregate Test Statistics for Global Null and associated pvalues
  avgstat0 = mean(gpstat0)
  avgstat0adj = mean(gpstat0adj)
  avgstat1 = mean(gpstat1)
  avgstat1adj = mean(gpstat1adj)

  pvavgstat0 = 1-stats::pnorm(avgstat0)
  pvavgstat0adj = 1-stats::pnorm(avgstat0adj)
  pvavgstat1 = 1-stats::pnorm(avgstat1)
  pvavgstat1adj = 1-stats::pnorm(avgstat1adj)

  #Key Player: Unconditional and Conditional on rejection of the global null

  #Key Player (unconditional)

  phatkp0 = which.max(gpstat0)
  phatkp1 = which.max(gpstat1)

  phatkp0adj = which.max(gpstat0adj)
  phatkp1adj = which.max(gpstat1adj)

  #Conditional

  if (pvavgstat0 < pvcutoffkp){
    phatkp0rej0 = phatkp0
  }else{
    phatkp0rej0 = 0
  }

  if (pvavgstat1 < pvcutoffkp){
    phatkp1rej0 = phatkp1
  }else{
    phatkp1rej0 = 0
  }

  if (pvavgstat0adj < pvcutoffkp){
    phatkp0adjrej0 = phatkp0adj
  }else{
    phatkp0adjrej0 = 0
  }

  if (pvavgstat1adj < pvcutoffkp){
    phatkp1adjrej0 = phatkp1adj
  }else{
    phatkp1adjrej0 = 0
  }

  return(list(avgstat0=avgstat0,pvavgstat0=pvavgstat0,
              avgstat1=avgstat1,pvavgstat1=pvavgstat1,
              avgstat0adj=avgstat0adj,pvavgstat0adj=pvavgstat0adj,
              avgstat1adj=avgstat1adj,pvavgstat1adj=pvavgstat1adj,
              phatkp0=phatkp0,phatkp1=phatkp1,
              phatkp0adj=phatkp0adj,phatkp1adj=phatkp1adj,
              phatkp0rej0=phatkp0rej0,phatkp1rej0=phatkp1rej0,
              phatkp0adjrej0=phatkp0adjrej0,phatkp1adjrej0=phatkp1adjrej0))

}
