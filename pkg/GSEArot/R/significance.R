significance <-
function(ES, ESnull, doplot=FALSE) { 
    
    #Include observed ES in null distribution to prevent p-values of magnitude 0
    ESnull <- c(ES,ESnull)
    #Calculate p-value separately for positive (>= 0) and negative ES as the 
    #number of ES's in the null distribution as least as extreme as the observed ES
    if(ES >= 0) {
        Nl <- sum(ESnull >= ES)
        pos <- ESnull[ESnull >= 0]
        p.value <- Nl/length(pos)
        NES <- ES/mean(pos)
        NESnull <- ESnull/mean(pos)
    } else {
        Ns <- sum(ESnull <= ES)
        neg <- ESnull[ESnull < 0]
        p.value <- Ns/length(neg)
        NES <- ES/abs(mean(neg))
        NESnull <- ESnull/abs(mean(neg))
    }
    #Remove observed NES from null distribution
    NESnull <- NESnull[-1]
    if(doplot)   
        hist(NESnull,freq=FALSE,xlab="NES",ylab="Density",main="Null distribution for NES",...)
    
    return(list(p.value=p.value, NES=NES, NESnull=NESnull))
}
