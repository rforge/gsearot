t.90th <-
function( X, y=NULL ) {
    
    if( is.null(y) ){
        #One-sample test
        n <- ncol(X)
        m <- apply(X,2,FUN=mean,na.rm=TRUE)
        std <- apply(X,2,FUN=sd,na.rm=TRUE)
        percStd <- quantile(std,probs=0.9)
        stat <- m / ( (std + percStd) / sqrt(2*n) )
    } else {   
        #Two-sample test 
        n1 <- sum(y==0)
        n2 <- sum(y)
        m1 <- apply(X[y==0,],2,FUN=mean,na.rm=TRUE)
        m2 <- apply(X[y==1,],2,FUN=mean,na.rm=TRUE)
        v1 <- apply(X[y==0,],2,FUN=var,na.rm=TRUE)
        v2 <- apply(X[y==1,],2,FUN=var,na.rm=TRUE)
        std <- sqrt(v1/n1 + v2/n2)
        percStd <- quantile(std,probs=0.9)
        stat <- (m1 - m2) / (std + percStd)
    }
        
    stat
}
