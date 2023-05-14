
#include "functions_and_data.stan"

parameters {
  real zlambda;
}

transformed parameters {
  vector[4] lpNBE[ngames];
  real<lower=0> lambda = exp(prior_lambda[1]+prior_lambda[2]*zlambda);
  for (gg in 1:ngames) {
      real pNU = exp(lpNash[gg][1]);
      real pNL = exp(lpNash[gg][3]);
      real b1 = UU1[gg][1,1]-UU1[gg][1,2]-UU1[gg][2,1]+UU1[gg][2,2];
      real b2 = UU2[gg][1,1]-UU2[gg][1,2]-UU2[gg][2,1]+UU2[gg][2,2];
      real pU = inv_logit(-signum(b1)*lambda*(pNU+signum(b2)*lambda*pNL)
      /(1-signum(b1*b2)*lambda^2));
      real pL = inv_logit(-signum(b2)*lambda*(pNL+signum(b1)*lambda*pNU)
      /(1-signum(b1*b2)*lambda^2));
              
      
      
      
      
      
      lpNBE[gg] = [log(pU),log(1-pU),log(pL),log(1-pL)]';
    }
  
  
}

model {
  
  
    for (cc in 1:ncohorts) {
    target += UseData[cc] * UDLRcount[cc]'*lpNBE[GameID[cc]];
  
    }
  zlambda ~ std_normal();
  
  
  
}

generated quantities {
  vector[ncohorts] log_lik;
  
    for (cc in 1:ncohorts) {
    log_lik[cc] = UDLRcount[cc]'*lpNBE[GameID[cc]];
  
  }
}

