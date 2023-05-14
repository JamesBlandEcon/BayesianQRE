
#include "functions_and_data.stan"

parameters {

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
      lpQRE[gg] = SolveLuceQRE(lambda*exp(taulambdaGame*zG[gg]),U,A,nc,ftol);
    }
    
    
  
  
  
}

model {
  
  zlambda ~ std_normal();
  ztaulambdaGame~ std_normal();
  zG ~ std_normal();
  
  for (cc in 1:ncohorts) {
    target += UseData[cc]*UDLRcount[cc]'*lpQRE[GameID[cc]];
  }
  
}

generated quantities {
  vector[ncohorts] log_lik;
  for (cc in 1:ncohorts) {
      log_lik[cc] = UDLRcount[cc]'*lpQRE[GameID[cc]];
    }
}


