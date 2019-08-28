jointA <- function(z, x, lambda1, lambda2, wc, alpha, k){
    if(z < (x / alpha)){
        return(0)
    } else{    
        tmp <- massA(x, lambda1, wc, alpha) * nonfeedingA(z-(x / alpha), lambda1, lambda2, wc, alpha, k)
        return(tmp)
    }
}

