reml_genedep <- function(alphaStart, Xlist, Ylist, Zlist, TsList, GsList, BsList, A, b) {
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
  #    - Bs: structure of batch dependencies between samples in a gene set (assuming same design for all gene sets).
  #    - set: vector of indices indicating gene set for each sample.
  #    - A and b are restrictions for parameters. See the constrOptim() help file.
  #Output:
  #    - $alpha.hat: vector of estimates for the parameters in alphaStart.
  #    - $beta.hat: vector of estimated parameters for fixed design factors.
  #    - $sigma: estimated common variance
  #    - $maklix: the maximum likelihood value
  
  p <- ncol(Xlist[[1]])             #Number of fixed factors
  nset <- length(Xlist)
  
  res <- constrOptim( theta=alphaStart, f=.loglik_genedep, grad=NULL, ui=A, ci=b, X=Xlist, Z=Zlist, Y=Ylist, Ts=TsList, Gs=GsList, Bs=BsList)
  alpha.hat <- res$par
  maxlik <- res$value
  
  phi.hat <- alpha.hat[3]
  ga.hat <- alpha.hat[4]
  #Time variance, gene variance, random error variance, remaining variances
  theta.hat <- c(1,alpha.hat[-c(3,4)])
  
  RSS <- numeric(nset) 
  N <- 0
  for(k in 1:nset) {
    
    Y <- Ylist[[k]]
    X <- Xlist[[k]]
    Z <- Zlist[[k]]
    Bs <- BsList[[k]]
    Gs <- GsList[[k]]
    Ts <- TsList[[k]]
    
    n <- length(Y)
    dim(Y) <- c(n,1)
    N <- N + n #Counting the total number of samples
    
    #Time and gene correlations
    Rt <- exp(-phi.hat*Ts)
    Rg <- exp(-ga.hat*Gs)
    
    #Make covariance matrix
    if( ncol(Z) > 2 ) {
      #only works when we have a batch effect. Else use makecov()
      V.hat <- theta.hat[4]*Bs + theta.hat[1]*Rt + theta.hat[2]*Rg + diag(n)*theta.hat[3]
    }else V.hat <- theta.hat[1]*Rt + theta.hat[2]*Rg + diag(n)*theta.hat[3]#V.hat <- makecov_genedep(theta.hat,Rt=Rt,Rg=Rg,design=Z)  
    
    Vinv <- solve(V.hat)
    XVinv <- t(X)%*%Vinv
    XVinvX <- XVinv%*%X
    XV <- solve(XVinvX)%*%XVinv
    
    beta.hat <- XV%*%Y
    r <- Y - X%*%beta.hat
    
    RSS[k] <- t(r)%*%Vinv%*%r        
  }    
  sigma.hat <- sum(RSS) / (N-nset*p)
  
  return(list(alpha=alpha.hat,beta=beta.hat,sigma=sigma.hat, maxlik=maxlik))
  
}