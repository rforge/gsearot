reml2 <-
function(alphaStart, X, Y, Z, Ts, Gs, Bs=NULL, set, A, b) {
  #Finds restricted maximum likelihood estimates of parameters.
  #
  #Input:
  #    - alphaStart: a vector of start values for the parameters to be estimated. The first parameter must be gene variance, 
  #      the second parameter random error variance, the third must be a parameter controlling correlation between time points and the fourth must
  #      be a parameter controlling correlation between genes. The remaining parameters can be variances for other random design factors. 
  #    - X: design matrix with fixed design variables for a gene set (assuming same design for all gene sets).
  #    - Y: vector of observed values for all gene sets.
  #    - Z: design matrix for random design factors for a gene set (assuming same design for all gene sets).
  #    - Ts: structure of time dependencies between samples in a gene set (assuming same design for all gene sets).
  #    - Gs: structure of gene dependencies between samples in a gene set (assuming same design for all gene sets).
  #    - Bs: (optional) structure of batch dependencies between samples in a gene set (assuming same design for all gene sets).
  #    - set: vector of indices indicating gene set for each sample.
  #    - A and b are restrictions for parameters. See the constrOptim() help file.
  #Output:
  #    - $alpha.hat: vector of estimates for the parameters in alphaStart.
  #    - $beta.hat: vector of estimated parameters for fixed design factors.
  #    - $sigma: estimated common variance
  #    - $maklix: the maximum likelihood value
  
  N <- length(Y)             #Number of samples
  p <- ncol(X)             #Number of fixed factors
  nset <- length(unique(set))
  n <- N/nset
  
  if(is.null(Bs)) Bs <- matrix(0,ncol(Ts),ncol(Gs))
  
  res <- constrOptim( theta=alphaStart, f=.loglik2, grad=NULL, ui=A, ci=b, X=X, Z=Z, Y=Y, Ts=Ts, Gs=Gs, Bs=Bs, set=set)
  alpha.hat <- res$par
  maxlik <- res$value
  
  phi <- alpha.hat[3]
  ga <- alpha.hat[4]
  #Time variance, gene variance, random error variance, remaining variances
  theta.hat <- c(1,alpha.hat[-c(3,4)])
  
  #Time and gene correlations
  Rt <- exp(-phi*Ts)
  Rg <- exp(-ga*Gs)
  
  #Make covariance matrix
  if( ncol(Z) > 2 ) {
    #only works when we have a batch effect. Else use makecov()
    V.hat <- theta.hat[4]*Bs + theta.hat[1]*Rt + theta.hat[2]*Rg + diag(n)*theta.hat[3]
  }else V.hat <- makecov2(theta.hat,Rt=Rt,Rg=Rg,design=Z) 
  
  Vinv <- solve(V.hat)
  XVinv <- t(X)%*%Vinv
  XVinvX <- XVinv%*%X
  XV <- solve(XVinvX)%*%XVinv
  
  RSS <- numeric(nset) 
  for(k in 1:nset) {
    
    idx <- which(set==k)
    Yk <- Y[idx]
    
    beta.hat <- XV%*%Yk
    r <- Yk - X%*%beta.hat
    RSS[k] <- t(r)%*%Vinv%*%r
  }
  
  sigma.hat <- sum(RSS) / (nset*(n-p))
  
  return(list(alpha=alpha.hat,beta=beta.hat,sigma=sigma.hat, maxlik=maxlik))
  
}
