
#include "functions_and_data.stan"

parameters {
  real zlambda;
  real zr;
  real ztaulambdaGame;
  real ztaurGame;
  vector[2] zG[ngames];
}

transformed parameters {
  vector[4] lpQRE[ngames];
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
      
      
      for (ro in 1:4) {
        for (co in 1:4) {
          if (UU[gg][ro,co]==0) {
            U[ro,co]=0;
          } else {
            U[ro,co]=exp(exp(theta[2])*log(UU[gg][ro,co]));
          }
          
        }
      }
      lpQRE[gg] = SolveQRE(exp(theta[1]),U,A,nc,ftol);
      
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
      target += UseData[cc]*UDLRcount[cc]'*lpQRE[GameID[cc]];
    }
  
  
}

generated quantities {
vector[ncohorts] log_lik;
  for (cc in 1:ncohorts) {
      log_lik[cc] = UDLRcount[cc]'*lpQRE[GameID[cc]];
    }
}
