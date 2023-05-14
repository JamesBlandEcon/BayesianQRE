
#include "functions_and_data.stan"

parameters {
  real<lower=0> kappa;
  real zr;
}

transformed parameters {
  
  
  real<lower=0> r;
  vector[4] lpQRE[ncohorts];
  r = exp(prior_r[1]+zr*prior_r[2]);

  
  
    
    
    for (cc in 1:ncohorts) {
      vector[2] theta;
      matrix[4,4] U;
      
      
      for (ro in 1:4) {
        for (co in 1:4) {
          if (UU[GameID[cc]][ro,co]==0) {
            U[ro,co]=0;
          } else {
            U[ro,co]=exp(r*log(UU[GameID[cc]][ro,co]));
          }
          
        }
      }
     {
        vector[2] ptmp;
        vector[2] pN;
        ptmp = U[1:2,3:4]\[1,1]';
        pN[1] = ptmp[1]/sum(ptmp);
        ptmp = U[3:4,1:2]\[1,1]';
        pN[2] = ptmp[1]/sum(ptmp);
       lpQRE[cc] =  [log(pN[1]),log(1.0-pN[1]),log(pN[2]),log(1.0-pN[2])]';
     }
      //SolveQRE(lambda,U,A,nc,ftol);
      
    }

    
  
  
  
}

model {
  
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
    real pU = exp(lpQRE[cc][1]);
    real pL = exp(lpQRE[cc][3]);
    
    log_lik[cc]= ( lbeta(kappa*pU+UDLRcount[cc][1],kappa*(1.0-pU)+UDLRcount[cc][2])
                -lbeta(kappa*pU,kappa*(1.0-pU)));
    
    log_lik[cc] += ( lbeta(kappa*pL+UDLRcount[cc][3],kappa*(1.0-pL)+UDLRcount[cc][4])
                -lbeta(kappa*pL,kappa*(1.0-pL)));
  
  }
}
