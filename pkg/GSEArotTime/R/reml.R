reml <-
function(alphaStart, X, Y, Z, subject, A, b) {
  #Finds restricted maximum likelihood estimates of parameters.
  #
  #Input:
  #    - alphaStart: a vector of start values for the parameters to be estimated. The first parameter must be random error variance, the last
  #             must be a parameter controlling correlation between time points. The remaining parameters can be varainces
  #             for other random design factors. 
  #    - X: design matrix with fixed design variables.
  #    - Y: matrix of observed values. Rows are samples and columns are variables.
  #    - Z: design matrix for random design factors.
  #    - subject: vector of indices indicating subject for each sample.
  #    - A and b are restrictions for parameters. See the constrOptim() help file.
  #Output:
  #    - $alpha.hat: vector of estimates for the parameters in alphaStart.
  #    - $beta.hat: vector of estimated parameters for fixed design factors.
  #    - $sigma: estimated common variance
  #    - $maklix: the maximum likelihood value
  
  G <- ncol(Y) #Number of genes
  N <- nrow(Y) #Number of samples
  p <- ncol(X) #Number of fixed factors
  
  res <- constrOptim( theta=alphaStart, f=.loglik, grad=NULL, ui=A, ci=b, X=X, Z=Z, Y=Y, subject=subject)
  alpha.hat <- res$par
  maxlik <- res$value
  
  phi <- alpha.hat[length(alpha.hat)]
  theta.hat <- c(1,alpha.hat[-length(alpha.hat)])
  
  #Make covariance matrix with estimated parameters
  #Find correlation between time points
  R <- matrix(0,N,N)
  for(i in 1:(N-1)) {
    for(j in (i+1):N) {
      if(subject[i]==subject[j]) R[i,j] <- exp( -phi*( Z[j,1] - Z[i,1] ) )
    }
  }
  R <- R + t(R)
  diag(R) <- 1
  #Find covariance matrix
  V.hat <- makecov(theta.hat,correlation=R,design=Z,subject=subject)
  
  #Estimate beta and sigma for final values of alpha
  Vinv <- solve(V.hat)
  XVinv <- t(X)%*%Vinv
  XVinvX <- XVinv%*%X
  beta.hat <- solve(XVinvX)%*%XVinv%*%Y
  r <- Y - X%*%beta.hat
  RSS <- sum( diag( t(r)%*%Vinv%*%(r) ) )
  sigma.hat <- RSS / ( G*(N-p) )
  
  return(list(alpha=alpha.hat,beta=beta.hat,sigma=sigma.hat, maxlik=maxlik))
  
}
