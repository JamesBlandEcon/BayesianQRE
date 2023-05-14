
// include the previously shown function and data blocks here

#include "functions.stan"

#include "data.stan"

parameters {
  // median lambda. 
  real<lower=0> lambda;
  
  // NEW PARAMETERS ADDED IN HERE:
  
  // standard deviatoin of lambda
  real<lower=0> lambda_SD;
  // lambda specific to each game
  vector<lower=0>[ngames] lambda_g;
}



model {
  /* Compute QRE given lambda_g, and store the (log) equilibrium probabilities
  in lpQRE
  */
  vector[4] lpQRE[ngames];
  
  for (gg in 1:ngames) { // loop over the 12 games
      // extract the payoffs of game gg
      vector[8] payoffs = to_vector(Payoffs[gg,]);
      // compute QRE
      lpQRE[gg] = logitQRE(lambda_g[gg],payoffs,nc); // <--- note lambda replaced with lambda_g[gg]
      
      
    }
  
  
  // now loop over each cohort's choice frequency
    for (cc in 1:ncohorts) {
      // increment the log-likelihood 
      target += sum(UDLRcount[cc].*lpQRE[GameID[cc]]);
  
    }
    
  // specify the prior for lambda
    lambda ~ lognormal(prior_lambda[1],prior_lambda[2]);
    
  // NEW HIERARCHICAL STRUCTURE SPECIFIED HERE
  // prior for lambda_SD
    lambda_SD ~ lognormal(prior_lambda_sd_game[1],prior_lambda_sd_game[2]);
  
  // link between population parameters (lambda,lambda_SD) and the game-specific lambda_g
    lambda_g ~ lognormal(log(lambda),lambda_SD);
}


