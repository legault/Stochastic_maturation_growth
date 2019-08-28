feedingA <- function(x, lambda1, wc, alpha){
    if(x < (wc / alpha)){
        return(0)
    } else{
        tmp <- lambda1 * exp(-lambda1 * (x - (wc / alpha)))
        return(tmp)
    }
}

