nonfeedingA <- function(x, lambda1, lambda2, wc, alpha, k){
    if(x <= 0){
        return(0)
    } else{    
        n <- k / alpha # number of terms to sum
        temp.stor <- dim(n)
        for(i in 1:n){
            temp.stor[i] <- exp(i * log(lambda2) + (i - 1) * log(x) - (lambda2 * x) + nbouts(i, lambda1, wc, alpha, log = TRUE) - lfactorial(i - 1))
        }
        return(sum(temp.stor))
    }
}
