jointB <- function(z, x, j, u1, lambdaY, lambdaX, alpha){
    dgamma(z, j, u1) * feedingBpdf(x, z * alpha, lambdaX / alpha, lambdaY / alpha)  
}

