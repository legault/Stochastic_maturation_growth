###### Code for simulating part A of the caterpillar development - the critical weight period

# Create nine parameter sets
param <- expand.grid(lambdaX = c(1, 2, 4),
                     lambdaY = c(1, 2, 4),
                     wc = c(6),
                     alpha = c(.1))
# Empty list for storing simulation results
sim.stor <- list()
# Simulate non-feeding and feeding during part A
## Number of simulations
sims <- 100000
## Set seed
set.seed(20180914)
## Main loop
##
## Warning: Takes ~45 minutes with sims = 100000
##
for(i in 1:nrow(param)){
    sim.stor[[i]] <- matrix(NA, nrow = sims, ncol = 5)
    for(j in 1:sims){
        lambdaX <- param[i, 1]
        lambdaY <- param[i, 2]
        wc <- param[i, 3]
        alpha <- param[i, 4]
        counter <- 0
        feed <- 0
        mass <- 0
        nfeed <- 0
        while(mass < wc){
            counter <- counter + 1
            nfeed <- nfeed + rexp(1, lambdaY)
            t.feed <- rexp(1, lambdaX)
            feed <- feed + t.feed
            mass <- mass + (t.feed * alpha)
        }
        sim.stor[[i]][j, 1] <- counter
        sim.stor[[i]][j, 2] <- feed
        sim.stor[[i]][j, 3] <- nfeed
        sim.stor[[i]][j, 4] <- mass
    }
    sim.stor[[i]][, 5] <- sim.stor[[i]][, 2] + sim.stor[[i]][, 3]
}

## # Uncomment to save simulation results
## write(sim.stor[[1]][, c(4,5)], file = "Data/jointA-01.csv", sep = ",")
## write(sim.stor[[2]][, c(4,5)], file = "Data/jointA-02.csv", sep = ",")
## write(sim.stor[[3]][, c(4,5)], file = "Data/jointA-03.csv", sep = ",")
## write(sim.stor[[4]][, c(4,5)], file = "Data/jointA-04.csv", sep = ",")
## write(sim.stor[[5]][, c(4,5)], file = "Data/jointA-05.csv", sep = ",")
## write(sim.stor[[6]][, c(4,5)], file = "Data/jointA-06.csv", sep = ",")
## write(sim.stor[[7]][, c(4,5)], file = "Data/jointA-07.csv", sep = ",")
## write(sim.stor[[8]][, c(4,5)], file = "Data/jointA-08.csv", sep = ",")
## write(sim.stor[[9]][, c(4,5)], file = "Data/jointA-09.csv", sep = ",")

# Check nbout probability mass function (pmf)
## Source function and display value
source("pdfs/nbouts.R")$value
## Estimate probabilty mass of bout number from simulations
bout.num1 <- table(sim.stor[[1]][, 1]) / sims
bout.num2 <- table(sim.stor[[2]][, 1]) / sims
bout.num3 <- table(sim.stor[[3]][, 1]) / sims

## Plot simulation results with pmf solutions
### General plotting parameters
par(mfrow = c(1, 3), mar = c(5, 5, 2, 1), lend = 2)
col1 <- rgb(.9, 0, 0, .8) # color 1 (simulations)
col2 <- rgb(0, 0, .9, 1) # color 2 (pdfs)

### Panel A
plot(NA, xlim = c(40, 280), ylim = c(0, .06), bty = "n", yaxt = "n", xaxt = "n",
     xlab = "Number of feeding bouts", ylab = "Probability mass", cex.lab = 1.6)
axis(1, at = seq(40, 280, 40), tcl = -.5)
axis(2, at = seq(0, .06, .01), las = 1, tcl = -.5)
mtext(side = 3, expression(paste("(a) ", lambda[X], "=1", sep = "")), adj = 0, cex = 1.2)
### Add simulation results
clip(39.3, 280.5, -1, .06)
for(i in 1:length(bout.num1)){
    segments(as.numeric(names(bout.num1[i])), 0,
             as.numeric(names(bout.num1[i])), bout.num1[[i]],
             lwd = 1, col = col1)
}
### Solve pmf and add to plot
sol <- nbouts(seq(40, 280, length.out = 80), param[1, 1], param[1, 3], param[1, 4])
points(sol ~ seq(40, 280, length.out = 80), type = "p", pch = 1, col = col2, cex = 1)
legend("topright", c("Simulations", expression(paste("pmf for ", italic("K")[A]))), pch = c(15, 1), col = c(col1, col2), bty = "n", inset = 0.07, cex = 1.6)

### Panel B
plot(NA, xlim = c(40, 280), ylim = c(0, .06), bty = "n", yaxt = "n", xaxt = "n",
     xlab = "Number of feeding bouts", ylab = "", cex.lab = 1.6)
axis(1, at = seq(40, 280, 40), tcl = -.5)
axis(2, at = seq(0, .06, .01), las = 1, tcl = -.5)
### Add simulation results
clip(39.3, 280.5, -1, .06)
for(i in 1:length(bout.num2)){
    segments(as.numeric(names(bout.num2[i])), 0,
             as.numeric(names(bout.num2[i])), bout.num2[[i]],
             lwd = 1, col = col1)
}
### Solve pmf and add to plot
sol <- nbouts(seq(40, 280, length.out = 80), param[2, 1], param[2, 3], param[2, 4])
points(sol ~ seq(40, 280, length.out = 80), type = "p", pch = 1, col = col2, cex = 1)
mtext(side = 3, expression(paste("(b) ", lambda[X], "=2", sep = "")), adj = 0, cex = 1.2)

### Panel C
plot(NA, xlim = c(40, 280), ylim = c(0, .06), bty = "n", yaxt = "n", xaxt = "n",
     xlab = "Number of feeding bouts", ylab = "", cex.lab = 1.6)
axis(1, at = seq(40, 280, 40), tcl = -.5)
axis(2, at = seq(0, .06, .01), las = 1, tcl = -.5)
### Add simulation results
clip(39.3, 280.5, -1, .06)
for(i in 1:length(bout.num3)){
    segments(as.numeric(names(bout.num3[i])), 0,
             as.numeric(names(bout.num3[i])), bout.num3[[i]],
             lwd = 1, col = col1)
}
### Solve pmf and add to plot
sol <- nbouts(seq(40, 280, length.out = 80), param[3, 1], param[3, 3], param[3, 4])
points(sol ~ seq(40, 280, length.out = 80), type = "p", pch = 1, col = col2, cex = 1)
mtext(side = 3, expression(paste("(c) ", lambda[X], "=4", sep = "")), adj = 0, cex = 1.2)


