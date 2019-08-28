# README

This repository contains R scripts associated with the paper:

*A stochastic model for predicting age and mass at maturity of insects*

The main purpose of these scripts is allow users to solve the joint distribution of the stochastic model described in the paper.

## simulate-partA.R

This script contains code for simulating mass gain and development time during the critical weight period. In addition, it compare simulation results to solutions of probability mass function for number of feeding bouts ("pdfs/nbout").

## simulate-partB.R

This script contains code for simulating mass gain and development time during the cessation of growth period. It uses Gillespie's stochastic simulation algorithm to model the degradation of juvenile hormone. In addition, it compare simulated development times with those expected by the model.

## solve-sensitivity.R

This script contains code that solves (numerically) the joint density of age and mass at maturity at three different values of lambda_X (1, 2, 4). The script relies heavily on the probability functions in folder "pdfs".


## pdfs/

This folder contains all probability mass / density functions associated with the development model. Users interested in solving for the final joint density (i.e., maturation and growth in periods A and B together) will want to use function "jointAB", which relies on all other functions in this folder.

## simu/

This folder contains functions associated with simulating the development model. Function "gillespie_full" is the primary function. Users interested in simulating the development model when parameter values are nonhomogeneous should consider swapping the inter-event sampler "hpp.R" with its nonhomogeneous analog "nhpp" found in the [following repository](https://github.com/legault/SSAplus) 
