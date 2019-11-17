feedingBpdf <- function(x, t, lambdaX, lambdaY){
    if(x <= 0 | x >= t){
        return(0)
    } else{    
    int1 <- integrate(f = function(H){
        ((besselI(2 * sqrt(lambdaX * lambdaY * H * (t - x)), 2) + besselI(2 * sqrt(lambdaX * lambdaY * H * (t - x)), 0)) * sqrt(H)) / (exp(lambdaX * H) * sqrt(lambdaX * lambdaY * H * (t - x)))}, lower = 0, upper = x)$value
    int2 <- integrate(f = function(H){
        besselI(2 * sqrt(lambdaX * lambdaY * H * (t - x)), 1) / (exp(lambdaX * H) * sqrt(H))}, lower = 0, upper = x)$value
    temp <- exp(-lambdaY * (t - x)) * (((besselI(2 * sqrt(lambdaX * lambdaY * (t - x) * x), 1) / (exp(lambdaX * x) * sqrt(x)) - (lambdaX * lambdaY * int1) / 2) * sqrt(lambdaX * lambdaY * (t - x))) - ((lambdaX * lambdaY * int2) / (2 * sqrt(lambdaX * lambdaY * (t - x))))) + (lambdaY * exp(-lambdaY * (t - x)) * ((int2) * sqrt(lambdaX * lambdaY * (t - x)) + 1))
    return(temp)
    }
}
