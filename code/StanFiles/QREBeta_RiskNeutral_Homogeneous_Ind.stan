
#include "functions_and_data.stan"

parameters {
  real zlambda;
  real<lower = 0> kappa;
}

transformed parameters {
  vector[4] lpQRE[ngames];
  real<lower=0> lambda = exp(prior_lambda[1]+prior_lambda[2]*zlambda);
 for (gg in 1:ngames) {
   matrix[4,4] U;
      for (ro in 1:4) {
        for (co in 1:4) {
          if (UU[gg][ro,co]==0) {
            U[ro,co]=0;
          } else {
            U[ro,co]=(UU[gg][ro,co]);
          }
          
        }
      }
      lpQRE[gg] = SolveQRE(lambda,U,A,nc,ftol);
    }
    
  
  
}

model {
  
  
  for (cc in 1:ncohorts) {
    real pU = exp(lpQRE[GameID[cc]][1]);
    real pL = exp(lpQRE[GameID[cc]][3]);
    
    target += UseData[cc]*( lbeta(kappa*pU+UDLRcount[cc][1],kappa*(1.0-pU)+UDLRcount[cc][2])
                -lbeta(kappa*pU,kappa*(1.0-pU)));
    
    target += UseData[cc]*( lbeta(kappa*pL+UDLRcount[cc][3],kappa*(1.0-pL)+UDLRcount[cc][4])
                -lbeta(kappa*pL,kappa*(1.0-pL)));
  
    }
  kappa ~ exponential(prior_kappa);  
  zlambda ~ std_normal();
  
  
  
  
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

