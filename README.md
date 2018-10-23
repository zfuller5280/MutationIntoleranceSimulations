# Mutation Intolerance Simulations
-----------------------------------
### This repository contains code used to run forward simulations under a constant population size model and a model of plausible demographic history for protein truncating variants (PTVs) from Fuller et al. (2018) "Measuring Intolerance to Mutation in Human Genetics"

### Simulation Details
The basis for the forward simulation code was originally developed in Simons et al. (2014) and further modified in Amorim et al. (2017) and Simons et al. (2018). The key difference here is that, instead of keeping track of frequencies of alleles, these simulations keep track of the distinct number of segregating sites in a finite number of mutational opportunities M. Each mutational opportunity is a biallelic site in a diploid individual, and we assume no intragenic recombination. At each site, mutations arise at a rate u and only arise on a background currently free of other PTV mutations. Each generation is formed by Wright-Fisher sampling with selection modeled by choosing parents according to their fitnesses, with a selection coefficient s and dominance coefficient h.

In the constant population size model, the number of generations is equal to 10N, where N is the effective population size. In the demographic model of population size changes in Europeans, the simulations have a burn-in period of 150k generations with a constant population size of 14,448. The first population size change occurs 55,940 generations ago. As in Simons et al. (2018), the population size changes and the generations they occur at are determined by piecing together the four haplotype MSMC for times corresponding to <170Kya and the two haplotype MSMC for more ancient times >170Kya inferred by Schiffels and Durbin (2016) from European (CEU) HapMap individuals. Finally, in the present generation the number of individuals is randomly sampled to match the number of non-Finnish Europeans (NFE) individuals in EXaC.

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

Note: to run large a number of replicates, it is recommended to parallelize the simulations to run on multiple processors or as array jobs on a cluster computing environment.


#### Methods for the Figure Legend from Fuller et al. (2018)

##### 1A
This panel shows the cumulative probability and density of pLI for two parameter combinations (blue line: s=0.1,h=0.9; red line: s=0.9,h=0.1) for a hypothetical gene with 225 PTV mutational opportunities and a mutation rate of 1.5e-8 under a constant population size of 100,000 individuals. To simulate these parameters:

````
PTV_count_simulations 1.5e-8 225 0.1 0.9 1000000 0 100000
PTV_count_simulations 1.5e-8 225 0.9 0.1 1000000 0 100000 
````
##### 1B
This panel is identical to 1B, except that instead of a constant population size, the demographic model of Schiffels & Durbin was used. To simulate these parameters:

````
PTV_count_simulations 1.5e-8 225 0.1 0.9 1000000 1 0
PTV_count_simulations 1.5e-8 225 0.9 0.1 1000000 1 0
````

##### 1C
This shows the probability of observing a PTV count of 3 (generated from a single simulation with s=0.1,h=0.9 and the same mutational and demographic parameters as 1B) for a grid of h and s values. To generate a single "observed" PTV count:

````
PTV_count_simulations 1.5e-8 225 0.1 0.9 1 1 0
````
Then, 1000000 simulations were run for each s and h parameter combination in the grid s=[0.01,0.02...0.99,1] and h=[0,0.01...0.99,1]. The output of each simulation was stored in a 3d Python numpy array, which was then used to estimate the likelihood of observing a PTV count of 3 across the parameter space using the script ptv_sim_lkhd.py.

##### 1D
This panel shows the behavior of pLI as a function of hs. The three lines correspond to three different gene lenths of PTV mutational opportunities (cyan=112, purple=225, yellow=550). The calculation of pLI requires an "expected" number of PTVs. To obtain a neutral expected number of PTVs for each gene length, we ran 1000000 simulation replicates with s=0,h=0:

````
PTV_count_simulations 1.5e-8 112 0 0 1000000 1 0
PTV_count_simulations 1.5e-8 225 0 0 1000000 1 0
PTV_count_simulations 1.5e-8 550 0 0 1000000 1 0
````
We took the mean number of PTVs from the set of neutral simulation replicates as the expected number of PTVs for each gene length. Then, similar to 1C, 5000000 simulations were run for each s and h parameter combination in the grid of values s=[0.01,0.02...0.99,1] and h=[0.01,0.02...0.99,1]. The simulation output for identical values of hs (e.g. s=0.1,h=.5 & s=0.5,h=.1) were concatenated together to bring the total number of simulations to 1000000. For each value of hs, we then calculated pLI for each "observed" replicate count of PTVs contained in the output using the R function 'calc_pli()'.

##### 1E
This panel depicts the behavior of pLI for a single value of hs, in this case s=0.1, h=.05 (hs=0.005). The blue histogram shows the distribution of PTV counts for these selection parameters for a gene with 225 PTV mutational opportunities, u=1.5e-8, and the demographic model of Schiffels & Durbin. To generate this distribution:
````
PTV_count_simulations 1.5e-8 225 0.1 0.5 1000000 1 0
````
To generate the expected number of PTVs under neutrality required by the calculation of pLI, a similar set of simulations were run, but with h=0 and s=0:
````
PTV_count_simulations 1.5e-8 225 0 0 1000000 1 0
````
The mean of this distribution (PTVs=18) was then used as the expected number of PTVs in the calculation of pLI. The red line shows what the value of pLI is for possible observed PTV counts for this expected value. To generate these values, the following R code was used:
````
pli_scores<-vector()
count=1
for(i in seq(0,20)){
  pli_scores[count]<-calc_pli(i,18)
  count=count+1
}
````
The inset in the plot shows the density of pLI scores calculated for each replicate in the simulations.
##### 1F
The final inset shows the distribution of PTV counts for three different selection parameter combinations of a gene with 225 PTV mutational opportunities, u=1.5e-8, and the Schiffels-Durbin demographic model. The selection paramters correspond to a neutral gene (s=0,h=0), completely recessive (s=.1,h=0), and a weakly selected dominant gene (s=.001,h=1). To generate the distributions for the three parameter combinations:

````
PTV_count_simulations 1.5e-8 225 0 0 1000000 1 0
PTV_count_simulations 1.5e-8 225 0.1 0 1000000 1 0
PTV_count_simulations 1.5e-8 225 0.001 1 1000000 1 0
````
