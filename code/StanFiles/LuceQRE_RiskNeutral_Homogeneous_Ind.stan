
#include "functions_and_data.stan"

parameters {
  real zlambda;
}

transformed parameters {
  vector[4] lpQRE[ngames];
  real<lower=0> lambda = exp(prior_lambda[1]+prior_lambda[2]*zlambda);
  for (gg in 1:ngames) {
      matrix[4,4] U=UU[gg];
      lpQRE[gg] = SolveLuceQRE(lambda,U,A,nc,ftol);
    }
  
  
}

model {
  
  
    for (cc in 1:ncohorts) {
    target += UseData[cc] * UDLRcount[cc]'*lpQRE[GameID[cc]];
  
    }
  zlambda ~ std_normal();
  
  
  
}

generated quantities {
  vector[ncohorts] log_lik;
  
    for (cc in 1:ncohorts) {
    log_lik[cc] = UDLRcount[cc]'*lpQRE[GameID[cc]];
  
  }
}

