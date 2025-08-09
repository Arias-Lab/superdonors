#!/usr/bin/env Rscript
# ------------------------------------------
# Date:        2025-08-05
# Description: This script is a Monte Carlo simulation tool requiring published datasets of alleles and corresponding frequency - probability 
# ------------------------------------------

#output files
#montecarlo_simulation.txt

#input files
#it needs as input 3 independent datasets containing the following structure
#column 1 are alleles
#column 2 are frequencies - probabilities
#dataset 1 is a dataframe containing hla allele frequencies
#dataset 2 is a dataframe containing hlb allele frequencies
#dataset 3 is a dataframe containing hlc allele frequencies
#the frequenices, as expected, must have a sum = 1

#three example datasets are included
#MCA_example_HLAA.txt
#MCA_example_HLAB.txt
#MCA_example_HLAC.txt

#example datasets
hlaa <- read.table(file = "MCA_example_HLAA.txt", header = TRUE)
hlab <- read.table(file = "MCA_example_HLAB.txt", header = TRUE)
hlac <- read.table(file = "MCA_example_HLAC.txt", header = TRUE)
#set number of individual haplotypes to simulate with Monte Carlo
x_monte_carlo_individuals <- 15
# set seed
seed <- set.seed(as.numeric(Sys.time()))
set.seed(seed)

#functions

# Sample n values from a category and probabilities pair
sample_with_probs <- function(values, probs, n) {
  s <- sum(probs)
  probs <- probs / s  # normalize probability defensively for literature loss of 1 by decimal loss in reporting
  sample(values, size = n, replace = TRUE, prob = probs)
}

simulate_hla <- function(n, hlaa, hlab, hlac) {
  output <- data.frame(
    sample = seq_len(n),
    HLAA1  = sample_with_probs(hlaa$HLA_A, hlaa$Frequency_Ahere, n),
    HLAA2  = sample_with_probs(hlaa$HLA_A, hlaa$Frequency_Ahere, n),
    HLAB1  = sample_with_probs(hlab$HLA_B, hlab$Frequency_Bhere, n),
    HLAB2  = sample_with_probs(hlab$HLA_B, hlab$Frequency_Bhere, n),
    HLAC1  = sample_with_probs(hlac$HLA_C, hlac$Frequency_Chere, n),
    HLAC2  = sample_with_probs(hlac$HLA_C, hlac$Frequency_Chere, n),
    stringsAsFactors = FALSE
  )
  output[] <- lapply(output, function(col) if (is.factor(col)) as.character(col) else col)
  output
}

#main

montecarlo_simulation <- simulate_hla(n=x_monte_carlo_individuals, hlaa = hlaa, hlab = hlab, hlac = hlac)
write.table(montecarlo_simulation, file = "montecarlo_simulation.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)

rm(hlaa, hlab, hlac, montecarlo_simulation, sample_with_probs, seed, simulate_hla, x_monte_carlo_individuals)
