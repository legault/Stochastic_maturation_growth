feedingBpdf <- function(x, t, lambda1, lambda2){
    if(x <= 0 | x >= t){
        return(0)
    } else{    
    int1 <- integrate(f = function(H){
        ((besselI(2 * sqrt(lambda2 * lambda1 * H * (t - (t - x))), 2) + besselI(2 * sqrt(lambda2 * lambda1 * H * (t - (t - x))), 0)) * sqrt(H)) / (exp(lambda2 * H) * sqrt(lambda2 * lambda1 * H * (t - (t - x))))}, lower = 0, upper = (t - x))$value
    int2 <- integrate(f = function(H){
        besselI(2 * sqrt(lambda2 * lambda1 * H * (t - (t - x))), 1) / (exp(lambda2 * H) * sqrt(H))}, lower = 0, upper = (t - x))$value
    temp <- exp(-lambda1 * (t - (t - x))) * (((besselI(2 * sqrt(lambda2 * lambda1 * (t - (t - x)) * (t - x)), 1) / (exp(lambda2 * (t - x)) * sqrt((t - x))) - (lambda2 * lambda1 * int1) / 2) * sqrt(lambda2 * lambda1 * (t - (t - x)))) - ((lambda2 * lambda1 * int2) / (2 * sqrt(lambda2 * lambda1 * (t - (t - x)))))) + (lambda1 * exp(-lambda1 * (t - (t - x))) * ((int2) * sqrt(lambda2 * lambda1 * (t - (t - x))) + 1))
    return(temp)
    }
}
