
#include "functions_and_data.stan"

parameters {
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

    
    target += UseData[cc]*lpQRE[cc]'*UDLRcount[cc];
  
    }

  
}

generated quantities {
vector[ncohorts] log_lik;
  for (cc in 1:ncohorts) {
    real pU = exp(lpQRE[GameID[cc]][1]);
    real pL = exp(lpQRE[GameID[cc]][3]);
    
    log_lik[cc]= lpQRE[cc]'*UDLRcount[cc];
  
  }
}
