jointAB <- function(z, x, lambdaX, lambdaY, wc, alpha, j, u1, k=50){
    integrate(Vectorize(function(H, z, x, lambdaX, lambdaY, wc, alpha, j, u1, k){
        nonfeedingA(z - (wc / alpha) - H, lambdaX, lambdaY, wc, alpha, k) * jointB(H, x - wc, j, u1, lambdaX, lambdaY, alpha)
    }),
    z = z,
    x = x,
    lambdaX = lambdaX,
    lambdaY = lambdaY,
    wc = wc,
    alpha = alpha,
    j = j,
    u1 = u1,
    k = k,
    lower = 0, upper = (z - (wc / alpha)))$value
}

