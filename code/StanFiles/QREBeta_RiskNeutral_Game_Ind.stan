
#include "functions_and_data.stan"

parameters {
  real<lower=0> kappa;
  vector[ngames] zG;
  real zlambda;
  real ztaulambdaGame;
}

transformed parameters {
  real<lower=0> lambda;
  real<lower=0> taulambdaGame;
  vector[4] lpQRE[ngames];
  lambda = exp(prior_lambda[1]+prior_lambda[2]*zlambda);
  taulambdaGame = exp(prior_taulambdaGame[1]+ztaulambdaGame*prior_taulambdaGame[2]);
  
  
    
    for (gg in 1:ngames) {
      matrix[4,4] U=UU[gg];
      lpQRE[gg] = SolveQRE(lambda*exp(taulambdaGame*zG[gg]),U,A,nc,ftol);
    }
    
    
  
  
  
}

model {
  
  zlambda ~ std_normal();
  ztaulambdaGame~ std_normal();
  zG ~ std_normal();
  
  for (cc in 1:ncohorts) {
    real pU = exp(lpQRE[GameID[cc]][1]);
    real pL = exp(lpQRE[GameID[cc]][3]);
    
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
    real pU = exp(lpQRE[GameID[cc]][1]);
    real pL = exp(lpQRE[GameID[cc]][3]);
    
    log_lik[cc]= ( lbeta(kappa*pU+UDLRcount[cc][1],kappa*(1.0-pU)+UDLRcount[cc][2])
                -lbeta(kappa*pU,kappa*(1.0-pU)));
    
    log_lik[cc] += ( lbeta(kappa*pL+UDLRcount[cc][3],kappa*(1.0-pL)+UDLRcount[cc][4])
                -lbeta(kappa*pL,kappa*(1.0-pL)));
  
  }
}


