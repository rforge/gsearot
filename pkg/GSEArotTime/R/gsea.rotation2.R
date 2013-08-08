gsea.rotation2 <-
function(S, Y, Xlist, Zlist, contrast, covmatList=NULL, nrot=10000, ES.p=1) {
  #GSEA rotation for longitudinal data
  #Assumes common covariance matrix for within gene sets
  
  #Input:
  #    - S: matrix of 0/1 indicating if a gene belongs to a gene set or not. 
  #      Each column represent a gene set, each row represent a gene
  #    - Y: gene expression data. A list of vectors, one for each gene set
  #    - Xlist: design data. A list of vectors, one for each gene set.
  #    - contrast: matrix of contrasts. Each column represent one contrast
  #    - covmatList: covariance between samples (optional). A list of vectors, one for each gene set
  #    - nrot: number of permutations (default: 10000)
  #    - Es.p=1: weighting parameter for enrichment score (default=1)
  #    
  #Output:
  #    - $Ngenes: size of gene sets
  #    - $ES: enrichment score
  #    - $NES: normalised enrichment score
  #    - $p.value: p.value per gene set
  #    - $q.value: q-value per gene set 
  
  require(limma)
  k <- ncol(contrast) #Number of interesting contrasts
  idx <- which( contrast==1, arr.ind=TRUE)[,1] #Get the index of interesting contrasts
  g <- nrow(S) #Number of genes
  n <- length(Y)/g #Number of samples per gene (assuming same number of samples)
  p <- nrow(contrast) #Number of design factors
  d <- n - p #Number of independent samples to calculate variance from
  p0 <- p-k #Number of uninteresting effects
  nset <- ncol(S)
  
  m <- 0
  l <- 0
  eff <- matrix(0,n,g)
  for(i in 1:nset) {     
    
    Xi <- Xlist[[i]]  
    Zi <- Zlist[[i]]
    yi <- Y[(m+1):(m+nrow(Xi))]
    dim(yi) <- c(nrow(Xi),1)
    m <- m + nrow(Xi)
    
    #Transform data to remove correlation between samples
    if( !is.null( covmatList ) ) {
      R <- chol(covmatList[[i]]) 
      Xi <- backsolve(R, Xi, transpose=TRUE)  
      yi <- backsolve(R, yi, transpose=TRUE) 
    }
    
    Xi <- cbind(Xi[,-idx],Xi[,idx]) #Arrange design matrix so interesting columns are last
    
    gi <- length(unique(Zi[,2])) #Number of genes in gene set
    
    for(j in 1:gi) {
      
      l <- l + 1
      geneid <- which(Zi[,2]==j) 
      Xj <- Xi[geneid,]
      yj <- yi[geneid]    
      
      QRdata <- qr(Xj)
      Q <- qr.Q(QRdata,complete=TRUE)  #Adding n-p random orthogonal columns
      eff[,l] <- t(Q)%*%yj  #Project y onto orthonormal basis for X. Projected onto the last n-p random orthogonal columns we achieve independent residuals.
    }
  } #nset
  
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
  dim(B) <- c(k,g) 
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
  modfr <- matrix(0,nrot,g)
  esrot <- matrix(0,nrot,nset)
  for(i in 1:nrot) {
    
    #Rotating y
    Z <- matrix(rnorm((d+k)^2),(d+k),(d+k))
    QRdata <- QR(Z)
    Qr <- QRdata$Q
    effr <- Qr %*% eff
    
    resr <- effr[(k+1):nrow(effr), , drop=FALSE] #n-p independent residuals to calculate variance from
    
    s2r <- colMeans(resr^2) #Estimate variance with independent residuals (last n-p columns)
    if (is.finite(d0)){ 
      sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d))
    } else{
      sdr.post <- sqrt(s02)
    }
    
    Br <- effr[1:k,] 
    dim(Br) <- c(k,g) 
    sdr.post.mat <- matrix(sdr.post,nrow(Br),ncol(Br),byrow=TRUE)
    modtr <- Br/sdr.post.mat#Modererte t-verdier.
    if( k > 1 ) {
      modfr[i,] <- apply(modtr,2,FUN=function(x) sum(x^2)/length(x)) #En moderert F-verdi for hvert gen (sum av kvadrerte t-verdier delt p? antall kontraster)d
    } else modfr[i,] <- modtr
    
    for (j in 1:nset) {
      esrot[i,j] <- es(modfr[i,],S[,j],p=ES.p)
    }
    
  }  
  
  #genewise p-values
  pvals <- numeric()
  for(j in 1:g) {
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
