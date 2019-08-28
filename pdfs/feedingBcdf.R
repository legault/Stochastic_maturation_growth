feedingBcdf <- function(x, t, lambda1, lambda2){
    exp(-lambda1 * (t - x)) * (1 + sqrt((lambda1 * lambda2) * (t - x)) * integrate(f = function(Y){
        exp(-lambda2 * Y) * (Y ^ -.5) * besselI(2 * sqrt((lambda1 * lambda2) * (t - x) * Y), 1, expon.scaled = TRUE) / exp(-2 * sqrt((lambda1 * lambda2) * (t - x) * Y))},
        lower = 0, upper = x)$value)
}
