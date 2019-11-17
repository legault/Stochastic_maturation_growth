# README

This repository contains R scripts associated with the paper:

*A stochastic model for predicting age and mass at maturity of insects*

The main purpose of these scripts is allow users to evaluate the joint distribution of the stochastic model described in the paper.

## simulate-evaluate.R

This script contains code for simulating mass gain and development time across both the critical weight period and cessation of growth period. For the cessation of growth period, JH degradation is simulated with GIllespie's stochastic algorithm (see folder /sim).

Following the simulations, there is code to evaluate the joint pdfs for age and mass at maturity.

Finally, contour plots are used to compare simulations with evaluations.

## pdfs/

This folder contains all probability mass / density functions associated with the development model. Users interested in solving for the final joint density (i.e., maturation and growth in periods A and B together) will want to use function "jointAB", which relies on all other functions in this folder.

## sim/

This folder contains functions associated with simulating the development model. Function "gillespie_full" is the primary function. Users interested in simulating the development model when parameter values are nonhomogeneous should consider swapping the inter-event sampler "hpp.R" with its nonhomogeneous analog "nhpp" found in the [following repository](https://github.com/legault/SSAplus) 
