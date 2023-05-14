
#include "functions_and_data.stan"

parameters {
  real zr;
  real ztaurGame;
  real zG[ngames];
}

transformed parameters {
  
  
  real<lower=0> r;
  real<lower=0> taurGame;
  vector[4] lpQRE[ncohorts];
  r = exp(prior_r[1]+zr*prior_r[2]);
  taurGame=exp(prior_taurGame[1]+ztaurGame*prior_taurGame[2]);
  
  
    
    
    for (cc in 1:ncohorts) {
      
      matrix[4,4] U;
      
      real rc = r*exp(taurGame*zG[GameID[cc]]);
      for (ro in 1:4) {
        for (co in 1:4) {
          if (UU[GameID[cc]][ro,co]==0) {
            U[ro,co]=0;
          } else {
            U[ro,co]=exp(rc*log(UU[GameID[cc]][ro,co]));
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
  ztaurGame ~ std_normal();
  zG ~ std_normal();

  for (cc in 1:ncohorts) {
    target += UseData[cc]*UDLRcount[cc]'*lpQRE[cc];
  }
  
}

generated quantities {
  vector[ncohorts] log_lik;
  for (cc in 1:ncohorts) {
      log_lik[cc] = UDLRcount[cc]'*lpQRE[cc];
    }
}
