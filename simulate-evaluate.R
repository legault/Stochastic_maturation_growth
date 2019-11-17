# Load libraries
library("MASS")
library("colorspace") # for good colors in joint density plots
library(parallel) # for parallel computing

# Make cluster
clus <- makeCluster(6) # adjust as appropriate

# Load parallel package into cluster
clusterEvalQ(clus, library(parallel))

# Source (and display) probability density functions
# I currently do not have a function for the development time as this is simply dgamma
source("pdfs/feedingBpdf.R")$value 
source("pdfs/jointB.R")$value
source("pdfs/jointAB.R")$value
source("pdfs/nonfeedingA.R")$value
source("pdfs/nbouts.R")$value

# Source functions for simulation stochastic compartment model
source("sim/gillespie_full.R")
source("sim/hpp.R")
source("sim/onecompartmentINT.R")
# Set parameters for compartment model
## Point process associated with degradation (see gillespie_full)
pproc <- matrix(c(-1), ncol = 1)
## Max time for simulation
maxtime <- 200

# Number of simulations
sims <- 1000
# Parameter values
param <- data.frame(lambdaX = c(1, 2, 4),
                    lambdaY = c(1),
                    wc = c(6),
                    alpha = c(.1),
                    j = c(100),
                    u1 = c(2))

# Simulate model for different parameters (not parallel)
#
## Set seed
set.seed(20180915)
## Empty list for storing simulation results
sim.stor <- list()
## Loop
for(i in 1:nrow(param)){
    print(i) # Current iteration of loop
    sim.stor[[i]] <- matrix(NA, nrow = sims, ncol = 6)
    for(j in 1:sims){
        lambdaX <- param[i, "lambdaX"]
        lambdaY <- param[i, "lambdaY"]
        wc <- param[i, "wc"]
        alpha <- param[i, "alpha"]
        counter <- 0
        # Part A
        feed <- 0
        mass <- 0
        nfeed <- 0
        while(mass < wc){
            counter <- counter + 1
            nfeed <- nfeed + rexp(1, lambdaY)
            t.feed <- rexp(1, lambdaX)
            feed <- feed + t.feed
            if(feed > wc / alpha){
                feed <- wc / alpha
                mass <- wc
            } else{
                mass <- mass + (t.feed * alpha)
            }
        }
        sim.stor[[i]][j, 1] <- nfeed + feed
        sim.stor[[i]][j, 2] <- mass
        # Part B
        temp <- gillespie_full(param[i, "j"], maxtime, list(u1 = param[i, "u1"]), onecompartmentINT, pproc, hpp)
        if(temp[nrow(temp), 2] == 0){
            sim.stor[[i]][j, 3] <- temp[nrow(temp), 1] # store exit time
        }
        else{
            sim.stor[[i]][j, 3] <- NA
        }
        tt <- 0 # restart switching clock
        maxtt <- sim.stor[[i]][j, 3]
        state <- 0
        feed <- 0
        while(tt < maxtt){
            if(state == 0){ # non-feeding
                temptt <- rexp(1, lambdaY)
                if((tt + temptt) > maxtt){
                    tt <- maxtt
                    break
                } else{               
                    tt <- tt + temptt
                    state <- 1
                }
            } else{ # feeding
                temptt <- rexp(1, lambdaX)
                if(tt + temptt > maxtt){
                    temptt <- (maxtt - tt)
                    feed <- feed + temptt
                    tt <- maxtt
                    break
                } else{
                    feed <- feed + temptt
                    tt <- tt + temptt
                    state <- 0
                }
            }
        }
        sim.stor[[i]][j, 4] <- feed * param[i, "alpha"]
    }
    sim.stor[[i]][, 5] <- sim.stor[[i]][, 1] + sim.stor[[i]][, 3]
    sim.stor[[i]][, 6] <- sim.stor[[i]][, 2] + sim.stor[[i]][, 4]
}

# Trim data, keeping only columns 5 and 6 (final development time and mass)
sim.stor[[1]] <- sim.stor[[1]][, c(5, 6)]
sim.stor[[2]] <- sim.stor[[2]][, c(5, 6)]
sim.stor[[3]] <- sim.stor[[3]][, c(5, 6)]

# Estimate joint density of the simulations
## Set grid size
gridsize <- 50 # adjust as needed
## Estimate density
### Parameter set 1
jointd1 <- kde2d(x = sim.stor[[1]][, 1],
                 y = sim.stor[[1]][, 2],
                 lims = c(140, 210, 7.5, 9.5),
                 n = gridsize)
### Parameter set 2
jointd2 <- kde2d(x = sim.stor[[2]][, 1],
                 y = sim.stor[[2]][, 2],
                 lims = c(200, 280, 7, 9),
                 n = gridsize)
### Parameter set 3
jointd3 <- kde2d(x = sim.stor[[3]][, 1],
                 y = sim.stor[[3]][, 2],
                 lims = c(300, 400, 6.5, 7.5),
                 n = gridsize)

# Evaluate the pdfs 
## WARNING: Takes ~20 minutes with only 1 core. For less time, consider a smaller grid size
## Export everything to cluster
clusterExport(clus, ls())
## Set 1
sol1 <- parSapplyLB(clus, 1:length(jointd1$x), function(i) { sapply(1:length(jointd1$y), function(j) {
    jointAB(z = jointd1$x[i],
            x = jointd1$y[j],
            lambdaX = param[1, "lambdaX"],
            lambdaY = param[1, "lambdaY"],
            wc = param[1, "wc"],
            alpha = param[1, "alpha"],
            j = param[1, "j"],
            u = param[1, "u1"])})})
### Transpose for plotting
sol1 <- t(sol1)

## Parameter Set 2
sol2 <- parSapplyLB(clus, 1:length(jointd2$x), function(i) { sapply(1:length(jointd2$y), function(j) {
    jointAB(z = jointd2$x[i],
            x = jointd2$y[j],
            lambdaX = param[2, "lambdaX"],
            lambdaY = param[2, "lambdaY"],
            wc = param[2, "wc"],
            alpha = param[2, "alpha"],
            j = param[2, "j"],
            u = param[2, "u1"])})})
### Transpose for plotting
sol2 <- t(sol2)

## Parameter Set 3
sol3 <- parSapplyLB(clus, 1:length(jointd3$x), function(i) { sapply(1:length(jointd3$y), function(j) {
    jointAB(z = jointd3$x[i],
            x = jointd3$y[j],
            lambdaX = param[3, "lambdaX"],
            lambdaY = param[3, "lambdaY"],
            wc = param[3, "wc"],
            alpha = param[3, "alpha"],
            j = param[3, "j"],
            u = param[3, "u1"])})})
### Transpose for plotting
sol3 <- t(sol3)

# Plots (recreating Figure A1 in the manuscript)
par(mfrow = c(3, 2), mar = c(5, 5, 2, 1), lend = 2)

# Panel A (Simulations, Parameter Set 1)
xmin <- 140
xmax <- 210
xdiff <- xmax - xmin
ymin <- 7.5
ymax <- 9.5
ydiff <- ymax - ymin
plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), bty = "n", yaxt = "n", xaxt = "n",
     xlab = expression(paste("Time (hours), ", italic("Z"))), ylab = expression(paste("Mass gain (g), ", italic("M"))))
axis(1, at = seq(xmin, xmax, 5), tcl = -.5)
axis(2, at = seq(ymin, ymax, .5), las = 1, tcl = -.5)
mtext(side = 3, expression(paste("(a) ", lambda[X], "=1, ", lambda[Y], "=1 (", bold("Simulations"), ")", sep = "")), adj = 0, cex = 1.2)
# Contours
clip(xmin + 0.05, xmax, ymin, ymax)
.filled.contour(jointd1$x, jointd1$y, jointd1$z, levels = seq(0, max(jointd1$z), length.out = 7), col = rev(heat_hcl(6, c = 0, l = c(0, 100), power = c(1/5, 1.3))))
rect(xmin + (xdiff * .8), head(seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), -1), xmin + (xdiff * .9), tail(seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), -1), col = rev(heat_hcl(7, c = 0, l = c(0, 100), power = c(1/5, 1.3))))
text(labels = format(round(seq(0, max(jointd1$z), length.out = 7), digits = 3), nsmall = 2), x = xmin + (xdiff * .95), y = seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), las = 2, cex = 1.1)
text(labels = c("Probability density"), x = xmin + (xdiff * .85), y = ymin + (ydiff * .975), cex = 1.1)

# Panel B (Evaluations, Parameter Set 1)
plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), bty = "n", yaxt = "n", xaxt = "n",
     xlab = expression(paste("Time (hours), ", italic("Z"))), ylab = "")
axis(1, at = seq(xmin, xmax, 5), tcl = -.5)
axis(2, at = seq(ymin, ymax, .5), las = 1, tcl = -.5)
mtext(side = 3, expression(paste("(b) ", lambda[X], "=1, ", lambda[Y], "=1 (", bold("pdf"), ")", sep = "")), adj = 0, cex = 1.2)
# Contours
clip(xmin + 0.05, xmax, ymin, ymax)
.filled.contour(jointd1$x, jointd1$y, sol1, levels = seq(0, max(sol1), length.out = 7), col = rev(heat_hcl(6, c = 0, l = c(0, 100), power = c(1/5, 1.3))))
rect(xmin + (xdiff * .8), head(seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), -1), xmin + (xdiff * .9), tail(seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), -1), col = rev(heat_hcl(7, c = 0, l = c(0, 100), power = c(1/5, 1.3))))
text(labels = format(round(seq(0, max(sol1), length.out = 7), digits = 3), nsmall = 2), x = xmin + (xdiff * .95), y = seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), las = 2, cex = 1.1)
text(labels = c("Probability density"), x = xmin + (xdiff * .85), y = ymin + (ydiff * .975), cex = 1.1)

# Panel C (Simulations, Parameter Set 2)
xmin <- 200
xmax <- 280
xdiff <- xmax - xmin
ymin <- 7
ymax <- 9
ydiff <- ymax - ymin
plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), bty = "n", yaxt = "n", xaxt = "n",
     xlab = expression(paste("Time (hours), ", italic("Z"))), ylab = expression(paste("Mass gain (g), ", italic("M"))))
axis(1, at = seq(xmin, xmax, 5), tcl = -.5)
axis(2, at = seq(ymin, ymax, .5), las = 1, tcl = -.5)
mtext(side = 3, expression(paste("(c) ", lambda[X], "=2, ", lambda[Y], "=1 (", bold("Simulations"), ")", sep = "")), adj = 0, cex = 1.2)
# Simulation contours
clip(xmin + 0.05, xmax, ymin, ymax)
.filled.contour(jointd2$x, jointd2$y, jointd2$z, levels = seq(0, max(jointd2$z), length.out = 7), col = rev(heat_hcl(6, c = 0, l = c(0, 100), power = c(1/5, 1.3))))
rect(xmin + (xdiff * .8), head(seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), -1), xmin + (xdiff * .9), tail(seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), -1), col = rev(heat_hcl(7, c = 0, l = c(0, 100), power = c(1/5, 1.3))))
text(labels = format(round(seq(0, max(jointd2$z), length.out = 7), digits = 3), nsmall = 2), x = xmin + (xdiff * .95), y = seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), las = 2, cex = 1.1)
text(labels = c("Probability density"), x = xmin + (xdiff * .85), y = ymin + (ydiff * .975), cex = 1.1)

# Panel D (Evaluations, Parameter Set 2)
plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), bty = "n", yaxt = "n", xaxt = "n",
     xlab = expression(paste("Time (hours), ", italic("Z"))), ylab = "")
axis(1, at = seq(xmin, xmax, 5), tcl = -.5)
axis(2, at = seq(ymin, ymax, .5), las = 1, tcl = -.5)
mtext(side = 3, expression(paste("(d) ", lambda[X], "=2, ", lambda[Y], "=1 (", bold("pdf"), ")", sep = "")), adj = 0, cex = 1.2)
# pdf contours
clip(xmin + 0.05, xmax, ymin, ymax)
.filled.contour(jointd2$x, jointd2$y, sol2, levels = seq(0, max(sol2), length.out = 7), col = rev(heat_hcl(6, c = 0, l = c(0, 100), power = c(1/5, 1.3))))
rect(xmin + (xdiff * .8), head(seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), -1), xmin + (xdiff * .9), tail(seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), -1), col = rev(heat_hcl(7, c = 0, l = c(0, 100), power = c(1/5, 1.3))))
text(labels = format(round(seq(0, max(sol2), length.out = 7), digits = 3), nsmall = 2), x = xmin + (xdiff * .95), y = seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), las = 2, cex = 1.1)
text(labels = c("Probability density"), x = xmin + (xdiff * .85), y = ymin + (ydiff * .975), cex = 1.1)

# Panel E (Simulations, Parameter Set 3)
xmin <- 300
xmax <- 400
xdiff <- xmax - xmin
ymin <- 6.5
ymax <- 8.5
ydiff <- ymax - ymin
plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), bty = "n", yaxt = "n", xaxt = "n",
     xlab = expression(paste("Time (hours), ", italic("Z"))), ylab = expression(paste("Mass gain (g), ", italic("M"))))
axis(1, at = seq(xmin, xmax, 5), tcl = -.5)
axis(2, at = seq(ymin, ymax, .5), las = 1, tcl = -.5)
mtext(side = 3, expression(paste("(e) ", lambda[X], "=4, ", lambda[Y], "=1 (", bold("Simulations"), ")", sep = "")), adj = 0, cex = 1.2)
# Simulation contours
clip(xmin + 0.05, xmax, ymin, ymax)
.filled.contour(jointd3$x, jointd3$y, jointd3$z, levels = seq(0, max(jointd3$z), length.out = 7), col = rev(heat_hcl(6, c = 0, l = c(0, 100), power = c(1/5, 1.3))))
rect(xmin + (xdiff * .8), head(seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), -1), xmin + (xdiff * .9), tail(seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), -1), col = rev(heat_hcl(7, c = 0, l = c(0, 100), power = c(1/5, 1.3))))
text(labels = format(round(seq(0, max(jointd3$z), length.out = 7), digits = 3), nsmall = 2), x = xmin + (xdiff * .95), y = seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), las = 2, cex = 1.1)
text(labels = c("Probability density"), x = xmin + (xdiff * .85), y = ymin + (ydiff * .975), cex = 1.1)

# Panel F (Evaluations, Parameter Set 3)
plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), bty = "n", yaxt = "n", xaxt = "n",
     xlab = expression(paste("Time (hours), ", italic("Z"))), ylab = "")
axis(1, at = seq(xmin, xmax, 5), tcl = -.5)
axis(2, at = seq(ymin, ymax, .5), las = 1, tcl = -.5)
mtext(side = 3, expression(paste("(f) ", lambda[X], "=4, ", lambda[Y], "=1 (", bold("pdf"), ")", sep = "")), adj = 0, cex = 1.2)
# pdf contours
clip(xmin + 0.05, xmax, ymin, ymax)
.filled.contour(jointd3$x, jointd3$y, sol3, levels = seq(0, max(sol3), length.out = 7), col = rev(heat_hcl(6, c = 0, l = c(0, 100), power = c(1/5, 1.3))))
rect(xmin + (xdiff * .8), head(seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), -1), xmin + (xdiff * .9), tail(seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), -1), col = rev(heat_hcl(7, c = 0, l = c(0, 100), power = c(1/5, 1.3))))
text(labels = format(round(seq(0, max(sol3), length.out = 7), digits = 3), nsmall = 2), x = xmin + (xdiff * .95), y = seq(ymin + (ydiff * .4), ymin + (ydiff * .9), length.out = 7), las = 2, cex = 1.1)
text(labels = c("Probability density"), x = xmin + (xdiff * .85), y = ymin + (ydiff * .975), cex = 1.1)
