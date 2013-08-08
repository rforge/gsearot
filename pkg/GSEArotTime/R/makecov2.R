makecov2 <-
function(theta, Rt, Rg, design) {
  #FOR BALANCED AND UNBALANCED DESIGNS.
  #Make covariance matrix for model with time series correlation, gene correlation and corelations due to other random effetcs.
  #
  #Input:
  #    - theta: a vector of variances for all random factors, where first element must be time series variance, second element must be gene variance,
  #      and third element must be the random error variance. The remaining elements can be variances for any other random design factors.
  #    - Rt: a matrix of time correlations between samples (usually zero for samples from different subjects).
  #    - Rg: a matrix of correlations between genes
  #    - design: nxp design matrix for the p random design factors. Must contain time in the first column and gene in the second. The time column must 
  #              be numeric and the time intervals must be of actual size.
  #Output:
  #    - V: covariance matrix.
  
  N <- nrow(design)
  
  V <- matrix(0,N,N)
  #Loop through each element of upper triangle
  for( i in 1:(N-1) ) {
    for( j in (i+1):N ) {
      V[i,j] <- theta[1]*Rt[i,j] + theta[2]*Rg[i,j]
      
      #Remaining factors in design matrix
      if( ncol(design) > 2 ) {
        for( k in 3:ncol(design) ) {
          if( design[i,k] == design[j,k] ) {
            V[i,j] <- V[i,j] + theta[k+1]
          }
        }
      }
    }
  }
  V <- V + t(V)
  diag(V) <- sum(theta) #All variances, including random error, appear on the diagonal
  
  V        
}
