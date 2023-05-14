// data.stan
data {
  int ncohorts; // number of cohorts
  int ngames; // number of games
  
  // Experiment data and parameters 
  matrix[ngames,8] Payoffs; // Payoffs
  vector[4] UDLRcount[ncohorts];
  int GameID[ncohorts]; // Integer identifiying each game
  
  
  real prior_lambda[2]; // prior for lambda (lognormal)
  real prior_lambda_sd_game[2]; // prior for standard deviation of lambda, only used in the heterogeneous model
  
  real prior_r[2]; // prior for CRRA parameter, only used in the CRRA model
  
  int nc; // number of corrector steps for computing QRE
  
  
}
