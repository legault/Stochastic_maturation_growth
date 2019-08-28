###### Code for simulating part B of the caterpillar development - the cessation of growth period

# Source functions for simulation of stochastic compartment model
source("sim/gillespie_full.R")
source("sim/hpp.R")
source("sim/onecompartmentINT.R")

# Set seed
set.seed(20180915)

## 9 parameter sets
param <- expand.grid(j = c(100), u1 = c(1, 2, 4), alpha = c(0.1))
## Point process associated with transition
pproc <- matrix(c(-1), ncol = 1)
## Max time
maxtime <- 200
## Run one simulation
gillespie_full(param[1, 1], maxtime, list(u1 = param[1, 2]), onecompartmentINT, pproc, hpp)

# Simulate 9 parameter sets
## Empty list for storage
sim.stor2 <- list()
## Number of simulations
sims <- 100000
## Simulation loop
##
## Warning: Will take > 45min with sims = 100000
##
for(i in 1:nrow(param)){
    sim.stor2[[i]] <- matrix(NA, nrow = sims, ncol = 1)
    for(j in 1:sims){
        temp <- gillespie_full(param[i, 1], maxtime, list(u1 = param[i, 2]), onecompartmentINT, pproc, hpp)
        if(temp[nrow(temp), 2] == 0){
            sim.stor2[[i]][j, 1] <- temp[nrow(temp), 1] # store exit time
        } else{
            sim.stor2[[i]][j, 1] <- NA
        }
    }
}

## # Uncomment to save simulation results
## write(sim.stor2[[1]][, ], file = "Data/devtimeB-1.csv", sep = ",")
## write(sim.stor2[[2]][, ], file = "Data/devtimeB-2.csv", sep = ",")
## write(sim.stor2[[3]][, ], file = "Data/devtimeB-3.csv", sep = ",")


# Check if simulation results for development time in part B line up with gamma distribution

## Plot simulation results with pdf solutions
### General plotting parameters
par(mfrow = c(1, 3), mar = c(5, 5, 2, 1), lend = 2)
col1 <- rgb(.9, 0, 0, .8) # color 1 (simulations)
col2 <- rgb(0, 0, .9, 1) # color 2 (pdfs)

### Panel A
plot(NA, xlim = c(0, 125), ylim = c(0, 0.16), bty = "n", yaxt = "n", xaxt = "n",
     xlab = "Time to degrade JH (hours)", ylab = "Probability density", cex.lab = 1.6)
axis(1, at = seq(0, 125, 25), tcl = -.5)
axis(2, at = seq(0, .16, .04), las = 1, tcl = -.5, pos = 0)
### Add simulation results
points(density(sim.stor2[[1]], from = 0, to = 125), type = "l", col = col1, lwd = 3)
### Solve pdf and add to plot
sequ <- seq(0, 125, length.out = 80)
sol <- dim(length(sequ))
for(j in 1:length(sequ)){
    sol[j] <- dgamma(sequ[j], shape = param[1, 1], rate = param[1, 2])
}
points(sol ~ sequ, type = "p", pch = 1, col = col2, cex = 2)
mtext(side = 3, expression(paste("(a) ", mu, "=1", sep = "")), adj = 0, cex = 1.2)
legend("topright", c("Simulations", expression(paste("pdf for ", italic("Z")[B]))), pch = c(NA, 1), lty = c("solid", NA), lwd = c(3, NA), col = c(col1, col2), bty = "n", inset = 0.07, cex = 1.6)

### Panel B
plot(NA, xlim = c(0, 125), ylim = c(0, .16), bty = "n", yaxt = "n", xaxt = "n",
     xlab = "Time to degrade JH (hours)", ylab = "", cex.lab = 1.6)
axis(1, at = seq(0, 125, 25), tcl = -.5)
axis(2, at = seq(0, 0.16, .04), las = 1, tcl = -.5, pos = 0)
### Add simulation results
points(density(sim.stor2[[2]], from = 0, to = 125), type = "l", col = col1, lwd = 3)
### Solve pdf and add to plot
sequ <- seq(0, 125, length.out = 80)
sol <- dim(length(sequ))
for(j in 1:length(sequ)){
    sol[j] <- dgamma(sequ[j], shape = param[2, 1], rate = param[2, 2])
}
points(sol ~ sequ, type = "p", pch = 1, col = col2, cex = 2)
mtext(side = 3, expression(paste("(b) ", mu, "=2", sep = "")), adj = 0, cex = 1.2)

### Panel C
plot(NA, xlim = c(0, 125), ylim = c(0, .16), bty = "n", yaxt = "n", xaxt = "n",
     xlab = "Time to degrade JH (hours)", ylab = "", cex.lab = 1.6, col = col1, lwd = 3)
axis(1, at = seq(0, 125, 25), tcl = -.5)
axis(2, at = seq(0, .16, .04), las = 1, tcl = -.5, pos = 0)
### Add simulation results
points(density(sim.stor2[[3]], from = 0, to = 125), type = "l", col = col1, lwd = 3)
### Solve pdf and add to plot
sequ <- seq(0, 125, length.out = 80)
sol <- dim(length(sequ))
for(j in 1:length(sequ)){
    sol[j] <- dgamma(sequ[j], shape = param[3, 1], rate = param[3, 2])
}
points(sol ~ sequ, type = "p", pch = 1, col = col2, cex = 2)
mtext(side = 3, expression(paste("(c) ", mu, "=4", sep = "")), adj = 0, cex = 1.2)
