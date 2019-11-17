nbouts <- function(k, lambdaX, wc, alpha, log = FALSE){
    tmp <- -lambdaX * (wc / alpha) + (k - 1) * log(lambdaX * (wc / alpha)) - lfactorial(k - 1)
    if(log == TRUE){
        return(tmp)
    } else{
        return(exp(tmp))
    }
}
