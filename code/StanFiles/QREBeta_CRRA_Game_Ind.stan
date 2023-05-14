
#include "functions_and_data.stan"

parameters {
  real<lower=0> kappa;
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
    real pU = exp(lpQRE[GameID[cc]][1]);
    real pL = exp(lpQRE[GameID[cc]][3]);
    
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
