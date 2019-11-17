## One-compartment death model
onecompartmentINT <- function(t, X, param){
    with(as.list(c(param)),{
        death1 <- ifelse(X[1] > 0, u1, 0)
        c(death1)})
}
