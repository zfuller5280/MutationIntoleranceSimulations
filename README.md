# Mutation Intolerance Simulations
-----------------------------------
### This repository contains code used to run forward simulations under a constant population size model and a model of plausible demographic history for protein truncating variants (PTVs) from Fuller et al. (2018) "Measuring Intolerance to Mutation in Human Genetics"

#### Requirements
Compiled with: 
 - Boost 1.67
 - gcc 4.9

#### To compile
```
g++ -O3 -std=c++11 PTV_count_simulations.cpp -o PTV_count_simulations
```

#### To run 

The simulation requires 7 parameters and must be in the following order: Mutation rate, number of PTV mutational opportunities, selection coefficient, dominance coefficient, number of simulation runs, option to run demographic model (0 for constant size, 1 for Schiffels-Durbin model), and effective population size (only applicable if a constant population size model is used, otherwise set to 0).

For example, to run 1000 simulations for a gene of typical length (225 PTV mutational opportunities) with a selection coefficient of 0.5 and dominance coefficient of 0.5 using the Schiffels-Durbin model for a mutation rate of 1.5e-8.
```
PTV_count_simulations 1.5e-8 225 0.5 0.5 1000 1 0
```
Or, to run 1000 simulations for a gene of typical length (225 PTV mutational opportunities) with a selection coefficient of 0.5 and dominance coefficient of 0.5 under a constant population size model of Ne=100,000 for a mutation rate of 1.5e-8.
```
PTV_count_simulations 1.5e-8 225 0.5 0.5 1000 0 100000
```
The simulation writes to stdout, and each line contains the number of segregating PTV counts in the population matching the number of Non-Finnish Europeans contained in EXaC (33370 individuals) from the end of each simulation run. 

#### Methods for the Figure Legend from Fuller et al. (2018)

##### 1A
This panel show the cumulative probability and density of pLI for two parameter combinations (blue line: s=0.1,h=0.9; red line: s=0.9,h=0.1) for a hypothetical gene with 225 PTV mutational opportunities and a mutation rate of 1.5e-8 under a constant population size of 100,000 individuals. 
