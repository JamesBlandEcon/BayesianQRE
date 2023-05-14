
#include "functions_and_data.stan"

parameters {
  real zlambda;
  real zr;
}

transformed parameters {
  
  
  real<lower=0> lambda;
  real<lower=0> r;
  vector[4] lpNBE[ncohorts];
  lambda = exp(prior_lambda[1]+zlambda*prior_lambda[2]);
  r = exp(prior_r[1]+zr*prior_r[2]);

  
  
    
    
    for (cc in 1:ncohorts) {
      vector[2] theta;
      matrix[4,4] U;
      
      
      theta = [log(lambda),log(r)]';
      
       {
      real rc = exp(theta[2]);
      real l = exp(theta[1]);
      real pNU = exp(lpNash[GameID[cc]][1]);
      real pNL = exp(lpNash[GameID[cc]][3]);
      real b1 = UU1[GameID[cc]][1,1]^rc-UU1[GameID[cc]][1,2]^rc-UU1[GameID[cc]][2,1]^rc+UU1[GameID[cc]][2,2]^rc;
      real b2 = UU2[GameID[cc]][1,1]^rc-UU2[GameID[cc]][1,2]^rc-UU2[GameID[cc]][2,1]^rc+UU2[GameID[cc]][2,2]^rc;
      real pU = inv_logit(-signum(b1)*l*(pNU+signum(b2)*l*pNL)
      /(1-signum(b1*b2)*l^2));
      real pL = inv_logit(-signum(b2)*l*(pNL+signum(b1)*l*pNU)
      /(1-signum(b1*b2)*l^2));
      lpNBE[cc] = [log(pU),log(1-pU),log(pL),log(1-pL)]';
      }
      
      
    }
    
    
    
  
  
  
}

model {
  
  zlambda ~ std_normal();
  zr ~ std_normal();

  for (cc in 1:ncohorts) {

    
    target += UseData[cc]*lpNBE[cc]'*UDLRcount[cc];
  
    }

  
}

generated quantities {
vector[ncohorts] log_lik;
  for (cc in 1:ncohorts) {
    real pU = exp(lpNBE[GameID[cc]][1]);
    real pL = exp(lpNBE[GameID[cc]][3]);
    
    log_lik[cc]= lpNBE[cc]'*UDLRcount[cc];
  
  }
}
