## Inter-arrival time sampling function (homogenous Poisson process)
hpp <- function(intentemp){
    rexp(1, sum(intentemp))
}
