functions {
  
/* SolveQRE takes the model and game's parameters, and computes the logit QRE 
using corrector steps. if variable QRE ==0, then it computes the Nash equilibrium 
instead

The function returns log probabilities (log(pU),log(pD),log(pU),log(pD))
*/
vector SolveQRE(real lambda, matrix U,data matrix A,data int nc,data real ftol) {
     vector[2] q = [0,0]';
     
     
     for (cc in 1:nc) {
        vector[2] H;
        vector[4] p;
        matrix[2,2] DLq;
        p=1.0./(1+exp(-[q[1],-q[1],q[2],-q[2]]'));
        H=(q-lambda*A*U*p);
        if (sqrt(H'*H)>ftol) {  
        DLq = [[p[1]*p[2],0],[0,p[3]*p[4]]];
        q=q-([[1,0],[0,1]]-lambda*A*U*(A')*DLq)\H;
        }
     }
     
     return [-log(1+exp(-q[1])),-log(1+exp(q[1])),-log(1+exp(-q[2])),-log(1+exp(q[2]))]';
     
   }
//----------------------------------------------

/* SolveQRE2 allows lambda to be different between the two players
*/
vector SolveQRE2(vector lambda, matrix U,data matrix A,data int nc,data real ftol) {
     vector[2] q = [0,0]';
     
     
     for (cc in 1:nc) {
        vector[2] H;
        vector[4] p;
        matrix[2,2] DLq;
        p=1.0./(1+exp(-[q[1],-q[1],q[2],-q[2]]'));
        H=(q-(A*U*p).*lambda);
        if (sqrt(H'*H)>ftol) {  
        DLq = [[p[1]*p[2],0],[0,p[3]*p[4]]];
        q=q-([[1,0],[0,1]]-[[lambda[1],0],[0,lambda[2]]].*(A*U*(A')*DLq))\H;
        }
     }
     
     return [-log(1+exp(-q[1])),-log(1+exp(q[1])),-log(1+exp(-q[2])),-log(1+exp(q[2]))]';
     
   }

//----------------------------------------------
   
   
   // QRE with Luce errors - only works for 2x2 games
   vector SolveLuceQRE(real lambda, matrix U,data matrix A,data int nc,data real ftol) {
     vector[2] q = [0,0]';
     real a1 = U[1,3];
     real b1 = U[1,4];
     real c1 = U[2,3];
     real d1 = U[2,4];
     
     real a2 = U[3,1];
     real b2 = U[3,2];
     real c2 = U[4,1];
     real d2 = U[4,2];
     
     for (cc in 1:nc) {
        vector[2] H;
        vector[4] p;
        vector[4] u;
        matrix[2,2] DH;
        p=1.0./(1+exp(-[q[1],-q[1],q[2],-q[2]]'));
        u = U*p;
        
        H = [q[1]-lambda*(log(u[1])-log(u[2]))
              ,
               q[2]-lambda*(log(u[3])-log(u[4]))
              ]';
              
        DH = [
          [1.0         , -lambda*p[3]*(1-p[3])*((a1-b1)/u[1]-(c1-d1)/u[2])],
               [-lambda*p[1]*(1-p[1])*((a2-b2)/u[3]-(c2-d2)/u[4]), 1.0]
               ];
        
        if (sqrt(H'*H)>ftol) {  
        
        q=q-DH\H;
        }
     }
     
     return [-log(1+exp(-q[1])),-log(1+exp(q[1])),-log(1+exp(-q[2])),-log(1+exp(q[2]))]';
     
   }

//----------------------------------------------
   
   
   // QRE with Luce errors 2 lambdas
   vector SolveLuceQRE2(vector lambda, matrix U,data matrix A,data int nc,data real ftol) {
     vector[2] q = [0,0]';
     real a1 = U[1,3];
     real b1 = U[1,4];
     real c1 = U[2,3];
     real d1 = U[2,4];
     
     real a2 = U[3,1];
     real b2 = U[3,2];
     real c2 = U[4,1];
     real d2 = U[4,2];
     
     for (cc in 1:nc) {
        vector[2] H;
        vector[4] p;
        vector[4] u;
        matrix[2,2] DH;
        p=1.0./(1+exp(-[q[1],-q[1],q[2],-q[2]]'));
        u = U*p;
        
        H = [q[1]-lambda[1]*(log(u[1])-log(u[2]))
              ,
               q[2]-lambda[2]*(log(u[3])-log(u[4]))
              ]';
              
        DH = [
          [1.0         , -lambda[1]*p[3]*(1-p[3])*((a1-b1)/u[1]-(c1-d1)/u[2])],
               [-lambda[2]*p[1]*(1-p[1])*((a2-b2)/u[3]-(c2-d2)/u[4]), 1.0]
               ];
        
        if (sqrt(H'*H)>ftol) {  
        
        q=q-DH\H;
        }
     }
     
     return [-log(1+exp(-q[1])),-log(1+exp(q[1])),-log(1+exp(-q[2])),-log(1+exp(q[2]))]';
     
   }
   
   
  
   // Function returning the sign of its argument - needed for NBE
   int signum(real x) {
    return x < 0 ? -1 : x > 0;
  }
              
}
data {
  int ncohorts; // number of cohorts
  int ngames; // number of games
  int nZ; // number of random standard normals for Monte Carlo integration in HQRE
  
  // Experiment data and parameters 
  matrix[ngames,8] Payoffs; // Payoffs
  vector[4] UDLRcount[ncohorts];
  int GameID[ncohorts]; // Integer identifiying each game
  vector[ncohorts] UseData; // Does this cohort contribute to the likelihood function? Use for cross-validation
  //vector[2] pNash[ngames]; // Nash equilibrium probabilities
  
  
  real prior_lambda[2]; // prior for lambda (lognormal)
  real prior_r[2]; // prior for CRRA parameter r (lognormal)
  real prior_taulambdaGame[2];
  real prior_taurGame[2];
  real prior_taulambdaCohort[2];
  real prior_taurCohort[2];
  real prior_taulambdaBoth[2];
  real prior_taurBoth[2];
  
  
  real prior_lambda_sd_both[2];
  real prior_lambda_sd_game[2];
  real prior_lambda_sd_cohort[2];
  real prior_r_sd_both[2];
  real prior_r_sd_game[2];
  real prior_r_sd_cohort[2];
  real prior_lkj;
  real prior_lnSD[2];
  real prior_kappa; // Exponential
  
  int nc;
  real ftol;
  
}
transformed data {
  matrix[4,4] UU[ngames];
  matrix[2,4] A;
  matrix[2,2] UU1[ngames];
  matrix[2,2] UU2[ngames];
  vector[4] lpNash[ngames];
  vector[nZ] ZR = to_vector(normal_rng(rep_vector(0.0,nZ),1.0));
  int sizeMU;
  
  A = [[1,-1,0,0],[0,0,1,-1]];
  
  for (gg in 1:ngames) {
    vector[2] ptmp;
    vector[2] pN;
    matrix[4,4] U;
    vector[8] pay = to_vector(Payoffs[gg,]);
    UU[gg]=[[0,0,pay[1],pay[2]],[0,0,pay[3],pay[4]],[pay[5],pay[6],0,0],[pay[7],pay[8],0,0]];
    
    UU1[gg] = [[pay[1],pay[2]],[pay[3],pay[4]]];
    UU2[gg] = [[pay[5],pay[6]],[pay[7],pay[8]]];
    
    
    
    
        ptmp = UU[gg][1:2,3:4]\[1,1]';
        pN[1] = ptmp[1]/sum(ptmp);
        ptmp = UU[gg][3:4,1:2]\[1,1]';
        pN[2] = ptmp[1]/sum(ptmp);
        lpNash[gg] =  [log(pN[1]),log(1.0-pN[1]),log(pN[2]),log(1.0-pN[2])]';
    
  }
  
    
  
}
