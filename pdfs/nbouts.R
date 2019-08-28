nbouts <- function(k, lambda, wc, alpha, log = FALSE){
    tmp <- -lambda * (wc / alpha) + (k - 1) * log(lambda * (wc / alpha)) - lfactorial(k - 1)
    if(log == TRUE){
        return(tmp)
    } else{
        return(exp(tmp))
    }
}
