
#include "functions_and_data.stan"

parameters {
  real zlambda;
  real zr;
  real ztaulambdaGame;
  real ztaurGame;
  vector[2] zG[ngames];
}

transformed parameters {
  vector[4] lpNBE[ngames];
  real<lower=0> lambda = exp(prior_lambda[1]+zlambda*prior_lambda[2]);
  real<lower=0> r = exp(prior_r[1]+zr*prior_r[2]);
  real<lower=0> taulambdaGame=exp(prior_taulambdaGame[1]+ztaulambdaGame*prior_taulambdaGame[2]);
  real<lower=0> taurGame=exp(prior_taurGame[1]+ztaurGame*prior_taurGame[2]);
  
  
    {
    
    for (gg in 1:ngames) {
      vector[2] theta;
      matrix[4,4] U;
      
      vector[2] tauGame = [taulambdaGame,taurGame]';
      
      
      theta = [log(lambda),log(r)]'+tauGame .* zG[gg];
       {
      real rc = exp(theta[2]);
      real l = exp(theta[1]);
      real pNU = exp(lpNash[gg][1]);
      real pNL = exp(lpNash[gg][3]);
      real b1 = UU1[gg][1,1]^rc-UU1[gg][1,2]^rc-UU1[gg][2,1]^rc+UU1[gg][2,2]^rc;
      real b2 = UU2[gg][1,1]^rc-UU2[gg][1,2]^rc-UU2[gg][2,1]^rc+UU2[gg][2,2]^rc;
      real pU = inv_logit(-signum(b1)*l*(pNU+signum(b2)*l*pNL)
      /(1-signum(b1*b2)*l^2));
      real pL = inv_logit(-signum(b2)*l*(pNL+signum(b1)*l*pNU)
      /(1-signum(b1*b2)*l^2));
      lpNBE[gg] = [log(pU),log(1-pU),log(pL),log(1-pL)]';
      }
      
      
    }
    
    
    
  
  
  
}
}

model {
  
  zlambda ~ std_normal();
  zr ~ std_normal();
  ztaulambdaGame~ std_normal();
  ztaurGame~ std_normal();
  for (gg in 1:ngames) {
    zG[gg]~std_normal();
  }
  
  
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
