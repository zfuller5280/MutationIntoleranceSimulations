# Mutation Intolerance Simulations
-----------------------------------
### This repository contains code used to run forward simulations under a constant population size model and a model of plausible demographic history for protein truncating variants (PTVs) from Fuller et al. (2018) "Measuring Intolerance to Mutation in Human Genetics"

### Simulation Details
The basis for the forward simulation code was originally developed in Simons et al. (2014) and further modified in Amorim et al. (2017) and Simons et al. (2018) (see https://github.com/sellalab/ForwardSimulator). The key difference here is that, instead of keeping track of frequencies of alleles, these simulations keep track of the distinct number of segregating sites in a finite number of mutational opportunities *M*. Each mutational opportunity is a biallelic site in a diploid individual, and we assume no intragenic recombination. At each site, mutations arise at a rate *u* and only arise on a background currently free of other PTV mutations. Each generation is formed by Wright-Fisher sampling with selection modeled by choosing parents according to their fitnesses, with a selection coefficient *s* and dominance coefficient *h*.

In the constant population size model, the number of generations is equal to 10*N*, where *N* is the effective population size. In the demographic model of population size changes in Europeans, the simulations have a burn-in period of 150k generations with a constant population size of 14,448. The first population size change occurs 55,940 generations ago. As in Simons et al. (2018), the population size changes and the generations they occur at are determined by piecing together the four haplotype MSMC for times corresponding to <170Kya and the two haplotype MSMC for more ancient times >170Kya inferred by Schiffels and Durbin (2016) from European (CEU) HapMap individuals. Finally, in the present generation the number of individuals is randomly sampled to match the number of non-Finnish Europeans (NFE) individuals in EXaC.

#### Requirements
Compiled with: 
 - Boost 1.67
 - gcc 4.9

#### To compile
```
g++ -O3 -std=c++11 PTV_count_simulations.cpp -o PTV_count_simulations
```

#### To run 

The simulation requires 7 parameters and must be in the following order: Mutation rate, number of PTV mutational opportunities *M*, selection coefficient *s*, dominance coefficient *h*, number of simulation runs, option to run demographic model (0 for constant size, 1 for Schiffels-Durbin model), and effective population size *N* (only applicable if a constant population size model is used, otherwise set to 0).

For example, to run 1000 simulations for a gene of typical length (225 PTV mutational opportunities) with a selection coefficient of 0.5 and dominance coefficient of 0.5 using the Schiffels-Durbin model for a mutation rate of 1.5e-8.

```
PTV_count_simulations 1.5e-8 225 0.5 0.5 1000 1 0
```
Or, to run 1000 simulations for a gene of typical length (225 PTV mutational opportunities) with a selection coefficient of 0.5 and dominance coefficient of 0.5 under a constant population size model of *N*=100,000 for a mutation rate of 1.5e-8.
```
PTV_count_simulations 1.5e-8 225 0.5 0.5 1000 0 100000
```
The simulation writes to stdout, and each line contains the number of segregating PTV counts in the population matching the number of Non-Finnish Europeans contained in EXaC (33370 individuals) from the end of each simulation run.

Note: to run large a number of replicates, it is recommended to parallelize the simulations to run on multiple processors or as array jobs on a cluster computing environment.

