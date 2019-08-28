jointB <- function(z, x, j, u1, lambdaY, lambdaX, alpha){
    dgamma(z, j, u1) * massB(x, z, lambdaY, lambdaX, alpha)  
}

