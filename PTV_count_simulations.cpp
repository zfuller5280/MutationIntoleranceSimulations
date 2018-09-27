//Written by Zach Fuller for Fuller et al. (2018).
//Constant size simulations originally written as R script by Jeremy Berg.
//Demographic model incorporated using code from Yuval Simons for Simons et al. (2018)
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <map>
#include <time.h>
#include <string>
#include <iostream>
  using std::cout;
  using std::endl;
#include <iomanip>
  using std::setprecision;
#include <cstdlib>
  using std::atoi;
  using std::atof;
#include <ctime>
  using std::time;
#include <vector>
  using std::vector;
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/poisson_distribution.hpp>
  using boost::poisson_distribution;
#include <boost/random/binomial_distribution.hpp>
  using boost::binomial_distribution;
#include <boost/random/variate_generator.hpp>
  using boost::variate_generator;

//Mersenne twister RNG from Boost
boost::mt19937 gent;

int gens,RUNS,length;
int popsize(int demographic_model, int gen, int i);
int popNe;
double sel,DOM,mutU;
int idx=0;
int site_class;
int exp_lof_num;
double mut_rate,expec_lof,obs_num;
int initgen;
int demographic_model;

//Population sizse changes from Schiffles and Durbin (2015) inferred by MSMC of Europeans
int Ne[56]= {14448,14068,14068,14464,14464,15208,15208,16256,16256,17618,17618,19347,19347,21534,21534,24236,24236,27367,27367,30416,30416,32060,32060,31284,29404,26686,23261,18990,16490,16490,12958,12958,9827,9827,7477,7477,5791,5791,4670,4670,3841,3841,3372,3372,3287,3359,3570,4095,4713,5661,7540,11375,14310,13292,14522,613285};
int T[56]={55940,51395,47457,43984,40877,38067,35501,33141,30956,28922,27018,25231,23545,21951,20439,19000,17628,16318,15063,13859,12702,11590,10517,9482,8483,7516,6580,5672,5520,5156,4817,4500,4203,3922,3656,3404,3165,2936,2718,2509,2308,2116,1930,1752,1579,1413,1252,1096,945,798,656,517,383,252,124,0};

using namespace std;

//Functions for boost random poisson and binomial variables
int boost_poi(double l);
int boost_binom(double n, double p);

int main(int argc, char *argv[]){
      //Get the options. Must appear in the order of mutation rate, number of mutational opportunities, selection, dominance, runs, option to run demographic model, Ne if constant population size model is run
      if (argc==8)
         {
           mutU=atof(argv[1]);
           length=atof(argv[2]);
           sel=atof(argv[3]);
           DOM=atof(argv[4]);
           RUNS=atof(argv[5]);
           //Input 1 if Schiffles-Durbin demographic model should be used. Else, constant size is used
           demographic_model=atof(argv[6]);
           //Only useful if constant size population model is used
           popNe=atof(argv[7]);
         }
      else{
        cout<<"Error: Not enough parameters enetered" << endl;
      }
    gent.seed(time(NULL));

    for(int run=0;run<RUNS;run++){
      //Set the burn in (here total simulation length is 150000 generations for Schiffles-Durbin demographic model)
      if (demographic_model==1){
        //Burn-in period, only used under Schiffels-Durbin model
        initgen=150000;
      }
      else{
        //Burn-in period, only used under constant population size model
        initgen=popNe*10;
      }
      double freqs [length] = {};
      double counts [length] = {};
      double mut_counts [length] = {};
      double mut_freqs [length] = {};
      double sel_freqs [length] = {};
      //The intitial frequency
      double init_freq = 1./(2*popsize(demographic_model,initgen,0));
      double freq_sum = 0.0;
      int n_seg = 0;

      std::fill_n(freqs, length, init_freq);

      for(int i=0;i<length;i++){
        counts[i] = freqs[i]*2*popsize(demographic_model,initgen,0);
      }

      int gen=initgen;
      idx=0;
      while (gen>0){
        double total_freq=0.;
        for(int i=0;i<length;i++){
          mut_counts[i] = counts[i] + boost_poi(2*popsize(demographic_model,gen,i)*mutU*(1-freqs[i]));
          mut_freqs[i] = (mut_counts[i])/(2*popsize(demographic_model,gen,i));
          total_freq += mut_freqs[i];
        }

        //Calculate the fitness
        double wbar = (pow(total_freq,2)*(1-sel) + ((2*total_freq)*(1-total_freq)*(1-DOM*sel)) + (pow((1-total_freq),2)));

        for(int i=0;i<length;i++){
          sel_freqs[i] = (mut_freqs[i]*(total_freq*(1-sel)) + mut_freqs[i]*(1-total_freq)*(1-DOM*sel))/wbar;
          counts[i] = boost_binom(2*popsize(demographic_model,gen,i), sel_freqs[i]);
          freqs[i] = counts[i]/(2*popsize(demographic_model,gen,i));
        }
        gen--;
      }
      for(int i=0;i<length;i++){
        if(freqs[i]!=0){
          obs_num = boost_binom(60706, freqs[i]);
          if (obs_num > 0){
            freq_sum += freqs[i];
            n_seg += 1;
          }
        }
    }
    cout<<n_seg<<"\t"<<freq_sum<<endl;
    }
}

int boost_poi(double l) //Regular Poisson random variate (Boost)

{
    if (l==0.0){
      return 0.0;
    }
    else{
      boost::random::poisson_distribution<> dist(l);
      return dist(gent);
    }
}

int boost_binom(double n, double p) //Regular binomial random variate (Boost)

{
      boost::random::binomial_distribution<> dist(n, p);
      return dist(gent);
}

int popsize(int demographic_model,int gen, int i) //Calculates the size of population at generation gen according to Schiffles & Durbin's model
{
  if (demographic_model==1){
    if (gen<=T[idx]) idx+=1;
    return Ne[idx];
  }
  else{
    return popNe;
  }
}
