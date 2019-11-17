nonfeedingA <- function(x, lambdaX, lambdaY, wc, alpha, k){
    if(x <= 0){
        return(0)
    } else{    
        n <- k / alpha # number of terms to sum
        temp.stor <- dim(n)
        for(i in 1:n){
            temp.stor[i] <- exp(i * log(lambdaY) + (i - 1) * log(x) - (lambdaY * x) + nbouts(i, lambdaX, wc, alpha, log = TRUE) - lfactorial(i - 1))
        }
        return(sum(temp.stor))
    }
}
