// functions.stan
functions {

vector logitQRE(real lambda, vector payoffs, data int nc) {
     
     
     real a1 = payoffs[1];
     real b1 = payoffs[2];
     real c1 = payoffs[3];
     real d1 = payoffs[4];
     
     real a2 = payoffs[5];
     real b2 = payoffs[6];
     real c2 = payoffs[7];
     real d2 = payoffs[8];
     
     vector[2] l = to_vector([0,0]); // initial guess
     
     for (cc in 1:nc) {
        vector[2] H;
        vector[4] p;
        matrix[2,2] DH;
        
        p = 1.0./(1.0+exp(-l)); // convert logit probabilities to probability levels
        
        // objective function
        H[1] = l[1]-lambda*((a1-b1-c1+d1)*p[2]+b1-d1);
        H[2] = l[2]-lambda*((a2-b2-c2+d2)*p[1]+b2-d2);
        
        // Jacobian of objective function
        DH[1,1] = 1;
        DH[2,2] = 1;
        DH[1,2] = -lambda*(a1-b1-c1+d1)*p[2]*(1.0-p[2]);
        DH[2,1] = -lambda*(a2-b2-c2+d2)*p[1]*(1.0-p[1]);
        
        // corrector step
        l = l-DH\H;
        
        
     }
     
     // return the log QRE probabilities
     return [-log(1+exp(-l[1])),-log(1+exp(l[1])),-log(1+exp(-l[2])),-log(1+exp(l[2]))]';
     
   }
   
  
}
