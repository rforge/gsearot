gsea.rotation <-
function(S, y, X, contrast, covmat=NULL, nrot=10000, ES.p=1) {
  #GSEA rotation for longitudinal  data
  #Assumes common covariance matrix for all genes
  
  #Input:
  #    - S: matrix of 0/1 indicating if a gene belongs to a gene set or not. 
  #      Each column represent a gene set, each row represent a gene
  #    - y: matrix of gene expression data. Genes are columns and samples are rows (??)
  #    - X: design matrix. Each row represent a sample (?)
  #    - contrast: matrix of contrasts. Each column represent one contrast
  #    - covmat: covariance matrix (optional). 
  #    - nrot: number of permutations (default: 10000)
  #    - Es.p: weighting parameter for enrichment score (default=1)
  #    
  #Output:
  #    - $Ngenes: size of gene sets
  #    - $ES: enrichment score
  #    - $NES: normalised enrichment score
  #    - $p.value: p.value per gene set
  #    - $q.value: q-value per gene set 
  
  #install the limma package from Bioconductor if necessary
  suppressMessages(suppressWarnings(a <- require(limma)))
  if(!a){
    source("http://bioconductor.org/biocLite.R")
    biocLite("limma",ask=FALSE)
    library(limma)
  }

  k <- ncol(contrast) #Number of interesting contrasts
  idx <- which( contrast==1, arr.ind=TRUE)[,1] #Get the index of interesting contrasts
  n <- nrow(X)
  p <- ncol(X)
  d <- n - p #Number of independent samples to calculate variance from
  p0 <- p-k #Number of uninteresting effects
  nset <- ncol(S)
  
  #Transform data to remove correlation between samples
  if( !is.null( covmat ) ) {
    R <- chol(covmat)              
    y <- t(backsolve(R, t(y), transpose = TRUE))
    X <- backsolve(R, X, transpose = TRUE)     
  }
  
  X <- cbind(X[,-idx],X[,idx]) #Arrange design matrix so interesting columns are last
  
  QRdata <- qr(X)
  Q <- qr.Q(QRdata,complete=TRUE)  #Adding n-p random orthogonal columns
  eff <- t(Q)%*%t(y)           #Project y onto orthonormal basis for X. Projected onto the last n-p random orthogonal columns we achieve independent residuals.
  
  if(p0>0) {
    eff <- eff[-(1:p0),] #Remove first p-k uninteresting rows
  }          
  res <- eff[(k+1):nrow(eff), , drop=FALSE] #n-p independent residuals to calculate variance from
  
  s2 <- colMeans(res^2) #Estimate variance with independent residuals (last n-p columns)
  sv <- squeezeVar(s2, df = d)
  d0 <- sv$df.prior
  s02 <- sv$var.prior
  sd.post <- sqrt(sv$var.post)
  
  B <- eff[1:k,] 
  dim(B) <- c(k,nrow(y)) 
  sd.post.mat <- matrix(sd.post,nrow(B),ncol(B),byrow=TRUE) 
  modt <- B/sd.post.mat #Modererte t-verdier.
  #Compute F-values if testing more than one contrast, else keep t-values
  if( k > 1 ){
    modf <- apply(modt,2,FUN=function(x) sum(x^2)/length(x)) #En moderert F-verdi for hvert gen (sum av kvadrerte t-verdier delt p? antall kontraster)
  } else modf <- modt
  
  
  escore <- numeric(nset)
  for (j in 1:nset) {
    escore[j] <- es(modf,S[,j],p=ES.p)
  }
  
  #Rotation test
  modfr <- matrix(0,nrot,nrow(y))
  esrot <- matrix(0,nrot,nset)
  
  for(i in 1:nrot) {
    
    #Rotating y
    Z <- matrix(rnorm((d+k)^2),(d+k),(d+k))
    QRdata <- QR(Z)
    Qr <- QRdata$Q
    effr <- Qr %*% eff
    
    #effr <- rotation(eff,method=2)
    
    resr <- effr[(k+1):nrow(effr), , drop=FALSE] #n-p independent residuals to calculate variance from
    
    s2r <- colMeans(resr^2) #Estimate variance with independent residuals (last n-p columns)
    if (is.finite(d0)){ 
      sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
    } else{
      sdr.post <- sqrt(s02)
    }
    
    Br <- effr[1:k,] 
    dim(Br) <- c(k,nrow(y)) 
    sdr.post.mat <- matrix(sdr.post,nrow(Br),ncol(Br),byrow=TRUE)
    modtr <- Br/sdr.post.mat
    if( k > 1 ) {
      modfr[i,] <- apply(modtr,2,FUN=function(x) sum(x^2)/length(x))
    } else modfr[i,] <- modtr
    
    
    for (j in 1:nset) {
      esrot[i,j] <- es(modfr[i,],S[,j],p=ES.p)
    }
    
  }  
  
  #genewise p-values
  pvals <- numeric()
  for(j in 1:nrow(y)) {
    pvals[j] <- sum(modfr[,j]>=modf[j])/nrot
  } 
  
  #Estimate significance
  p.value <- numeric(nset)
  NES <- numeric(nset)
  NES.null <- matrix(0,nrot,nset)
  for( j in 1:nset) {
    sig <- significance(escore[j], esrot[,j])
    p.value[j] <- sig$p.value            #p-value
    NES[j] <- sig$NES                    #Normalised enrichment score
    NES.null[,j] <- sig$NESnull          #Normalised null distribution
  }
  
  
  #q-values
  #If the test statistic is F, use a "one-tail test", else use a "two-tail test"
  q.value <- numeric(nset)
  if(k > 1) {
    for(j in 1:nset) {
      q.value[j] <- ( ( sum(NES.null[,j] >= NES[j]) + 1)/( nrot + 1 ) ) / ( sum(NES >= NES[j])/nset ) 
    }
    q.value <- ifelse( q.value > 1, 1, q.value)
  } else q.value <- sapply(NES, FUN=fdr, NESobs=NES, NESnull=NES.null)
  
  
  len.set <- colSums(S)
  return(list(NGenes = len.set, ES=escore, NES=NES, p.value=p.value, q.value=q.value)) 
}
