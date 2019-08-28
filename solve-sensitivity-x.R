# Load parallel package
library(parallel)
# Make cluster
clus <- makeCluster(4) # adjust as appropriate
# Load parallel package into cluster
clusterEvalQ(clus, library(parallel))

# Part A pdfs
source("pdfs/nbouts.R")$value
source("pdfs/feedingA.R")$value
source("pdfs/nonfeedingA.R")$value
source("pdfs/massA.R")$value # Mass gain (part A)
source("pdfs/jointA.R")$value
# Part B pdfs
source("pdfs/feedingBcdf.R")$value
source("pdfs/derivB.R")$value
source("pdfs/massB.R")$value # pdf for mass (wrapper for derivB)
source("pdfs/jointB.R")$value
# Joint pdf for A and B combined (double integral)
source("pdfs/jointAB.R")$value

# Set parameters
## Effect of changing lambdaX; all other parameters constant
param <- data.frame(lambdaX = c(1, 2, 4),
                    lambdaY = c(2),
                    wc = c(6),
                    alpha = c(0.1),
                    j = c(100),
                    u1 = c(2))

# Z and M values to solve for
## Adjust length.out for larger grid size (10x10 currently)
Zvalues <- seq(85, 285, length.out = 10) 
Mvalues <- seq(6.9, 12.3, length.out = 10)

# Export everything to cluster
clusterExport(clus, ls())

# Run everything
#
# Warning: Solving will take > 2 hours for a 100x100 grid
#
## Set 1
sol1 <- parSapplyLB(clus, 1:length(Zvalues), function(i) { sapply(1:length(Mvalues), function(j) {
    jointAB(Zvalues[i], Mvalues[j], param[1, 1], param[1, 2], param[1, 3], param[1, 4], param[1, 5], param[1, 6])})})
## ## Uncomment if you wish to save
## write(sol1, file = "Data/results-x1.csv", sep = ",")

## Set 2
sol2 <- parSapplyLB(clus, 1:length(Zvalues), function(i) { sapply(1:length(Mvalues), function(j) {
    jointAB(Zvalues[i], Mvalues[j], param[2, 1], param[2, 2], param[2, 3], param[2, 4], param[2, 5], param[2, 6])})})
## ## Uncomment if you wish to save
## write(sol2, file = "Data/results-x2.csv", sep = ",")

## Set 3
sol3 <- parSapplyLB(clus, 1:length(Zvalues), function(i) { sapply(1:length(Mvalues), function(j) {
    jointAB(Zvalues[i], Mvalues[j], param[3, 1], param[3, 2], param[3, 3], param[3, 4], param[3, 5], param[3, 6])})})
## ## Uncomment if you wish to save
## write(sol3, file = "Data/results-x4.csv", sep = ",")
