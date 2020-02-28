jointB <- function(z, x, j, u1, lambda1, lambda2, alpha){
    dgamma(z, j, u1) * feedingBpdf(x, z * alpha, lambda1 / alpha, lambda2 / alpha)  
}

