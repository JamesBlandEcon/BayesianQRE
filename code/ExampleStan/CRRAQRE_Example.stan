
// include the previously shown function and data blocks here

#include "functions.stan"

#include "data.stan"

parameters {
  // our model has one parameter: lambda. We need to constrin it to be positive
  real<lower=0> lambda;
  
  real<lower=0> r; // <-- Declared new CRRA parameter here
}



model {
  /* Compute QRE given lambda, and store the (log) equilibrium probabilities
  in lpQRE
  */
  vector[4] lpQRE[ngames];
  
  for (gg in 1:ngames) { // loop over the 12 games
      // extract the payoffs of game gg
      vector[8] payoffs = to_vector(pow(Payoffs[gg,],r));
      // compute QRE
      lpQRE[gg] = logitQRE(lambda,payoffs,nc);
    }
  
  
  // now loop over each cohort's choice frequency
    for (cc in 1:ncohorts) {
      // increment the log-likelihood 
      target += sum(UDLRcount[cc].*lpQRE[GameID[cc]]);
  
    }
    
  // specify the prior for lambda
    lambda ~ lognormal(prior_lambda[1],prior_lambda[2]);
    
  // specify prior for r
    r ~ lognormal(prior_r[1],prior_r[2]);
  
  
  
}


