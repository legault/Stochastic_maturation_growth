jointAB <- function(z, x, lambda1, lambda2, wc, alpha, j, u1, k=50){
    inner <- function(z, x, lambda1, lambda2, wc, alpha, j, u1, k, I){
        integrate(Vectorize(function(H, z, x, lambda1, lambda2, wc, alpha, j, u1, k, I){
            jointA(z - H, I, lambda1, lambda2, wc, alpha, k) * jointB(H, x - I, j, u1, lambda2, lambda1, alpha)}),
            z = z,
            x = x,
            lambda1 = lambda1,
            lambda2 = lambda2,
            wc = wc,
            alpha = alpha,
            j = j,
            u1 = u1,
            k = k,
            I = I,
            lower = 0, upper = z - (wc / alpha))$value
    }
    integrate(Vectorize(function(z, x, lambda1, lambda2, wc, alpha, j, u1, k, I){
        inner(z, x, lambda1, lambda2, wc, alpha, j, u1, k, I)}),
        z = z,
        x = x,
        lambda1 = lambda1,
        lambda2 = lambda2,
        wc = wc,
        alpha = alpha,
        j = j,
        u1 = u1,
        k = k,
        lower = wc, upper = x)$value
}

