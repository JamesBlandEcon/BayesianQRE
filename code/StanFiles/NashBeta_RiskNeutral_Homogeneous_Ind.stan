
#include "functions_and_data.stan"

parameters {
  real<lower=0> kappa;
}

transformed parameters {

  
  
    

    
  
  
  
}

model {
  

  for (cc in 1:ncohorts) {
    real pU = exp(lpNash[GameID[cc]][1]);
    real pL = exp(lpNash[GameID[cc]][3]);
    
    target += UseData[cc]*( lbeta(kappa*pU+UDLRcount[cc][1],kappa*(1.0-pU)+UDLRcount[cc][2])
                -lbeta(kappa*pU,kappa*(1.0-pU)));
    
    target += UseData[cc]*( lbeta(kappa*pL+UDLRcount[cc][3],kappa*(1.0-pL)+UDLRcount[cc][4])
                -lbeta(kappa*pL,kappa*(1.0-pL)));
  
    }
  kappa ~ exponential(prior_kappa);  
  
}

generated quantities {
vector[ncohorts] log_lik;
   for (cc in 1:ncohorts) {
    real pU = exp(lpNash[GameID[cc]][1]);
    real pL = exp(lpNash[GameID[cc]][3]);
    
    log_lik[cc]= ( lbeta(kappa*pU+UDLRcount[cc][1],kappa*(1.0-pU)+UDLRcount[cc][2])
                -lbeta(kappa*pU,kappa*(1.0-pU)));
    
    log_lik[cc] += ( lbeta(kappa*pL+UDLRcount[cc][3],kappa*(1.0-pL)+UDLRcount[cc][4])
                -lbeta(kappa*pL,kappa*(1.0-pL)));
  
  }
}
