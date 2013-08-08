fdr <-
function(NES,NESobs,NESnull) { 
    nrot <- nrow(NESnull)
    
    Nnp <- numeric(nrot)
    Nnl <- numeric(nrot)
    Nnn <- numeric(nrot)
    Nns <- numeric(nrot)

    #For each permutation, counting the number of positive/negative NES and the number of NES more 
    #extreme than the given NES in the null distribution
    for(i in 1:nrot) {
               
        if(NES >= 0) {
            Nnp[i] <- sum(NESnull[i,] >= 0) + 1
            Nnl[i] <- sum(NESnull[i,] >= NES) + 1
        } else {
           Nnn[i] <- sum(NESnull[i,] < 0) + 1
           Nns[i] <- sum(NESnull[i,] <= NES) + 1
        }
    }
    
    #Counting the number of positive/negative NES and NES more extreme than the given NES,
    #among the NES of all gene sets to be tested. Finally, calculating the q-value.
    if(NES >= 0) {
    
        Np <- sum(NESobs >= 0)
        Nl <- sum(NESobs >= NES)
        Nnpl <- Nnl/Nnp 
        q <- mean(Nnpl)/(Nl/Np)
        
    } else {
        
        Nn <- sum(NESobs < 0)
        Ns <- sum(NESobs <= NES)
        Nnns <- Nns/Nnn
        q <- mean(Nnns)/(Ns/Nn)
        
    }
    if( q > 1) q <- 1
    q
}
