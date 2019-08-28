massB <- function(X, maxt, lambdaY, lambdaX, alpha){
    derivB(X, maxt * alpha, lambdaY / alpha, lambdaX / alpha)
}
