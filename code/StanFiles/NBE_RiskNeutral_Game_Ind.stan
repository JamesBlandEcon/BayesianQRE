
#include "functions_and_data.stan"

parameters {

  vector[ngames] zG;
  real zlambda;
  real ztaulambdaGame;
}

transformed parameters {
  real<lower=0> lambda;
  real<lower=0> taulambdaGame;
  vector[4] lpNBE[ngames];
  lambda = exp(prior_lambda[1]+prior_lambda[2]*zlambda);
  taulambdaGame = exp(prior_taulambdaGame[1]+ztaulambdaGame*prior_taulambdaGame[2]);
  
  
    
    for (gg in 1:ngames) {
       real l = lambda*exp(taulambdaGame*zG[gg]);
       real pNU = exp(lpNash[gg][1]);
      real pNL = exp(lpNash[gg][3]);
      real b1 = UU1[gg][1,1]-UU1[gg][1,2]-UU1[gg][2,1]+UU1[gg][2,2];
      real b2 = UU2[gg][1,1]-UU2[gg][1,2]-UU2[gg][2,1]+UU2[gg][2,2];
      real pU = inv_logit(-signum(b1)*l*(pNU+signum(b2)*l*pNL)
      /(1-signum(b1*b2)*l^2));
      real pL = inv_logit(-signum(b2)*l*(pNL+signum(b1)*l*pNU)
      /(1-signum(b1*b2)*l^2));
      lpNBE[gg] = [log(pU),log(1-pU),log(pL),log(1-pL)]';
    }
    
    
  
  
  
}

model {
  
  zlambda ~ std_normal();
  ztaulambdaGame~ std_normal();
  zG ~ std_normal();
  
  for (cc in 1:ncohorts) {
    target += UseData[cc]*UDLRcount[cc]'*lpNBE[GameID[cc]];
  }
  
}

generated quantities {
  vector[ncohorts] log_lik;
  for (cc in 1:ncohorts) {
      log_lik[cc] = UDLRcount[cc]'*lpNBE[GameID[cc]];
    }
}


