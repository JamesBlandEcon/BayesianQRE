
#include "functions_and_data.stan"

parameters {
  real<lower=0> kappa;
  real zlambda;
  real zr;
}

transformed parameters {
  
  
  real<lower=0> lambda;
  real<lower=0> r;
  vector[4] lpQRE[ncohorts];
  lambda = exp(prior_lambda[1]+zlambda*prior_lambda[2]);
  r = exp(prior_r[1]+zr*prior_r[2]);

  
  
    
    
    for (cc in 1:ncohorts) {
      vector[2] theta;
      matrix[4,4] U;
      
      
      theta = [log(lambda),log(r)]';
      
      
      for (ro in 1:4) {
        for (co in 1:4) {
          if (UU[GameID[cc]][ro,co]==0) {
            U[ro,co]=0;
          } else {
            U[ro,co]=exp(r*log(UU[GameID[cc]][ro,co]));
          }
          
        }
      }
      lpQRE[cc] = SolveLuceQRE(lambda,U,A,nc,ftol);
      
    }
    
    
    
  
  
  
}

model {
  
  zlambda ~ std_normal();
  zr ~ std_normal();

  for (cc in 1:ncohorts) {
    real pU = exp(lpQRE[cc][1]);
    real pL = exp(lpQRE[cc][3]);
    
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
