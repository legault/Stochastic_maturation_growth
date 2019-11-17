# SSA function (full accounting of events and times)
## Function `gillespie_full` takes as input the following (most defined above):
## * `init` = Array containing all initial conditions (e.g., initial population size)
## * `maxtime` = Desired end time of simulation
## * `param` = Array containing all intrinsic demographic parameters, rate functions, and the environment function required for calculating intensities in `inten`(below)
## * `inten` = Function which returns all intensities/probabilities of the point processes
## * `pproc` = Array containing the state changes caused by the point process (in the same order as `inten`; one row per state change)
## * `store_size` = The initial size of the storage array
gillespie_full <- function(init, maxtime, param, inten, pproc, hpp, store_size=10000) {
    N <- init
    nvars <- length(N)
    # Create storage matrix
    results <- matrix(nrow = store_size, ncol = nvars + 1 )
    colnames(results) <- c("t", paste(rep("N", nvars), 1:nvars, sep = ""))
    # Initialize
    tottime <- 0
    # Row index to keep track of rows in storage matrix        
    row <- 1
    while(tottime <= maxtime){
        # Expand storage if necessary
        if(row > nrow(results)){
            results <- rbind(results, matrix(NA, nrow = store_size, ncol = nvars + 1))
        }
        results[row,] <- c(tottime, N) # Store result
        intentemp <- inten(tottime, N, param) # Calculate intensity
        if(all(intentemp == 0)){
            break
        } else if(min(intentemp) < 0) {
            warning("Exiting with intensity less than 0")
            break
        } else {
            tau <- hpp(intentemp)
            tottime <- tottime + tau
            which.proc <- sample(1:nrow(pproc),
                          size = 1,
                          prob = intentemp) # Which process
            N <- N + pproc[which.proc, ] # Carry out events
            row <- row + 1
        }
    }
    if(tottime > maxtime){
        results[1:(row - 1), ]
    } else {
        results[1:row, ]
    }
}
