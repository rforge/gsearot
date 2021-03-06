rotation <-
function(Y, method) { 
    n <- dim(Y)[1]
    #Generate matrix with random standard normally distributed numbers
    W <- matrix(rnorm(n^2), nrow=n, ncol=n)
    
    if(method==1) {
        QRdata <- QR(W)
        Qm <- QRdata$Q
        
    }   else if(method==2) {
        E<-scale(W,scale=F)
        QRdata <- QR(E)
        W <- QRdata$Q[,1:(n-1)]
        E <- matrix(rnorm((n-1)^2),(n-1))
        QRdata <- QR(E)
        R.star <- QRdata$Q
        Qm <- n^(-1)*matrix(1,n,n)+(W%*%R.star)%*%t(W)
        
    }   else {
        stop("method must be either 1 or 2 \n")    
    }
        
    #Multiply rotation matrix with original data matrix
    X <- Qm%*%Y
    X
}
