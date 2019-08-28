massA <- function(x, lambda1, wc, alpha){
    tmp <- feedingA(x, lambda1 / alpha, wc, 1)
    return(tmp)
}
