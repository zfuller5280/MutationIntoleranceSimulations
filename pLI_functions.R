#Functions to calculate pLI, pREC, and pNULL from observed and expected PTV counts. The initial pi's used in the function are those obtained from Lek et al. (2016)
calc_pli<-function(obs, exp){
  init_pis<-c(0.2075797,0.4887326,0.3036877)
  lambdas<-c(1,0.463,0.089)
  expec_vec<-lambdas*exp
  pois_vec<-dpois(obs, lambda=expec_vec)
  init_prob<-(pois_vec * init_pis)/sum((pois_vec * init_pis))
  return(init_prob[3])
}

calc_prec<-function(obs, exp){
  init_pis<-c(0.2075797,0.4887326,0.3036877)
  lambdas<-c(1,0.463,0.089)
  expec_vec<-lambdas*exp
  pois_vec<-dpois(obs, lambda=expec_vec)
  init_prob<-(pois_vec * init_pis)/sum((pois_vec * init_pis))
  return(init_prob[2])
}

calc_pnull<-function(obs, exp){
  init_pis<-c(0.2075797,0.4887326,0.3036877)
  lambdas<-c(1,0.463,0.089)
  expec_vec<-lambdas*exp
  pois_vec<-dpois(obs, lambda=expec_vec)
  init_prob<-(pois_vec * init_pis)/sum((pois_vec * init_pis))
  return(init_prob[1])
}

