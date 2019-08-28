# Function for numerically solving the derivative of the feeding time cdf
derivB <- function(X, maxt, lambdaY, lambdaX){
    if(isTRUE(all.equal(X, 0))){
          dexp(maxt, lambdaY)
    } else{
        if(isTRUE(all.equal(X, maxt)) | X > maxt){
            0
        } else{
            x <- as.numeric(X)
            result <- numericDeriv(quote(feedingBcdf(x, maxt, lambdaY, lambdaX)), "x")
            attr(result, "gradient")[1]
        }
    }
}

