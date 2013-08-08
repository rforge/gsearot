gsea.rotation <-
function(X,y=NULL,S,nrot=500,fun=t.90th) {
    
    nset <- ncol(S)
    
    stat <- fun(X,y)       
        
    #Calculate an enrichment score for each gene set
    escore <- numeric(nset)
    for( i in 1:nset ) {
        escore[i] <- es(stat,S[,i])
    }
    
    #Calculate a null distribution for each gene set
    es.null <- matrix(nrow=nrot,ncol=nset)
    for(j in 1:nrot) {
    
        if( j%%100 == 0 ) cat(paste(j,"rotations completed... \n"))
        
       #Rotation method depends on the data type
        if( is.null(y) ) {
            Xrot <- rotation(X, method=1) 
        } else Xrot <- rotation(X, method=2)
                
        #Compute test statistics and enrichment scores for rotated data
        testRot <- fun(Xrot,y) 
        for(i in 1:nset) {   
            es.null[j,i] <- es(testRot,S[,i])
        }
    }
    
    #Estimating p-values and normalised enrichment scores
    pvals <- numeric(nset)
    NES <- numeric(nset)
    NES.null <- matrix(nrow=nrot,ncol=nset)
    for( i in 1:nset) {
        res <- significance(escore[i], es.null[,i])
        pvals[i] <- res$p.value            
        NES[i] <- res$NES                  
        NES.null[,i] <- res$NESnull        
    }
    
    #Estimating FDR q-values
    qvals <- sapply(NES, FUN=fdr, NESobs=NES, NESnull=NES.null)
    
    return(list(p.values=pvals,q.values=qvals))
}
