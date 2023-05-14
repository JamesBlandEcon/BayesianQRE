set.seed(42)
library(tidyverse)
library(rstan)
  # declare some options for Stan
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
d<-readRDS("SeltenChmura2008data.rds")
  
lambda<-1
SimSize<-1000

DecisionsPerGame<-c(50,400)

logitQRE<-function(lambda, # Choice precision
                   payoffs, # payoffs of the game
                   nc,# Number of corrector steps
                   ftol # function tolerance for inverting DH
){ 
  
  
  a1<-payoffs[1]
  b1<-payoffs[2]
  c1<-payoffs[3]
  d1<-payoffs[4]
  
  a2<-payoffs[5]
  b2<-payoffs[6]
  c2<-payoffs[7]
  d2<-payoffs[8]
  
  l<-c(0,0) # initial guess
  
  for (cc in 1:nc) {
    p<-1/(1+exp(-l)) # convert logit probabilities to probability levels
    
    # objective function
    H<-c(0,0)
    H[1]<-l[1]-lambda*((a1-b1-c1+d1)*p[2]+b1-d1)
    H[2]<-l[2]-lambda*((a2-b2-c2+d2)*p[1]+b2-d2)
    
    DH<-diag(c(1,1))
    DH[1,2]<- -lambda*(a1-b1-c1+d1)*p[2]*(1-p[2])
    DH[2,1]<- -lambda*(a2-b2-c2+d2)*p[1]*(1-p[1])
    
    # corrector step
    dl<--solve(DH, tol = ftol)%*%H
    l<-l+dl
    
    
  }
  
  # return the log QRE probabilities
  c(-log(1+exp(-l[1])),-log(1+exp(l[1])), -log(1+exp(-l[2])), -log(1+exp(l[2])))
  
}

Model<-stan_model(paste0("ExampleStan/QRE_Example.stan"),
                  auto_write=TRUE,
                  allow_undefined=TRUE,
                  verbose=FALSE,
                  includes=c("\n"))

CRRAModel<-stan_model(paste0("ExampleStan/CRRAQRE_Example.stan"),
                  auto_write=TRUE,
                  allow_undefined=TRUE,
                  verbose=FALSE,
                  includes=c("\n"))

PAYOFFS<-d$Payoffs

# Choice probs (pU,pL)
ChoiceProbs<-matrix(0,dim(PAYOFFS)[1],2)
for (gg in 1:d$ngames) {
  ChoiceProbs[gg,]<-exp(logitQRE(lambda,PAYOFFS[gg,],10,1e-17))[c(1,3)]
}
ChoiceProbs %>% print()

SimResults<-tibble()

for (dd in 1:length(DecisionsPerGame)) {
for (ss in 1:SimSize) {
  
  print(paste(dd,',',ss))
  
  decisions<-DecisionsPerGame[dd]
  
  # Draw from the distribution implied by logit QRE 
  UDLRcount<-matrix(0,dim(PAYOFFS)[1],4)
  for (gg in 1:d$ngames) {
    nU<-rbinom(1,decisions,ChoiceProbs[gg,1])
    nL<-rbinom(1,decisions,ChoiceProbs[gg,2])
    UDLRcount[gg,]<-c(nU,decisions-nU,nL,decisions-nL)
  }
  
  # update the data list
  
  d$UDLRcount<-UDLRcount
  d$GameID<-1:12
  d$ncohorts<-12
  
  #Estimate Bayesian model
  start_time<-Sys.time()
  BayesFit<-sampling(
    Model,
    data  = d,
    seed = 1234
  )
  end_time<-Sys.time()
  BayesFitTime<-(end_time-start_time)
  
  # Estimate CRRA Bayesian model
  start_time<-Sys.time()
  CRRABayesFit<-sampling(
    CRRAModel,
    data  = d,
    seed = 1234
  )
  end_time<-Sys.time()
  CRRABayesFitTime<-(end_time-start_time)
  
  # likelihood function for risk-neutral model
  neglogLike<-function(lambda) {
    
    #compute QRE given lambda, and store the (log) equilibrium probabilities
    #in lpQRE
    lpQRE<-c()
    for (gg in 1:d$ngames) {
      payoffs<-d$Payoffs[gg,]
      lpQRE<-rbind(lpQRE,logitQRE(lambda,payoffs,10,1e-17))
    }
    logL<-0
    #now loop over each cohort's choice frequency
    for (cc in 1:d$ncohorts) {
      logL<-logL+sum(d$UDLRcount[cc,]*lpQRE[d$GameID[cc],])
    }
    -logL
  }
  start_time<-Sys.time()
  MLEFit<-stats4::mle(minusl=neglogLike,start=list(lambda=1),
                      lower=list(lambda=0)
                      ,method="L-BFGS-B")
  end_time<-Sys.time()
  MLEFitTime<-(end_time-start_time)
  
  
  # likelihood function for CRRA model
  neglogLike<-function(lambda,r) {
    
    #compute QRE given lambda, and store the (log) equilibrium probabilities
    #in lpQRE
    lpQRE<-c()
    for (gg in 1:d$ngames) {
      payoffs<-d$Payoffs[gg,]^r
      lpQRE<-rbind(lpQRE,logitQRE(lambda,payoffs,10,1e-17))
    }
    logL<-0
    #now loop over each cohort's choice frequency
    for (cc in 1:d$ncohorts) {
      logL<-logL+sum(d$UDLRcount[cc,]*lpQRE[d$GameID[cc],])
    }
    -logL
  }
  start_time<-Sys.time()
  CRRAMLEFit<-stats4::mle(minusl=neglogLike,start=list(lambda=1,r=1),
                      lower=list(lambda=0,r=0)
                      ,method="L-BFGS-B")
  end_time<-Sys.time()
  CRRAMLEFitTime<-(end_time-start_time)
  
sumBayes<-summary(BayesFit)$summary
sumCRRABayes<-summary(CRRABayesFit)$summary

addThis<-tibble(
  lambda=lambda,
  decisions=decisions,
  lambda_MLE=MLEFit@fullcoef,
  lambda_MLEsd=MLEFit@vcov %>% sqrt() %>% as.numeric(),
  lambda_B=sumBayes["lambda","mean"],
  lambda_Bsd = sumBayes["lambda","sd"],
  r_MLE=-1,
  r_MLEsd=-1,
  r_B=-1,
  r_Bsd = -1,
  Btime = BayesFitTime,
  MLEtime = MLEFitTime
)

SimResults<-rbind(SimResults,
                  addThis)

addThis<-tibble(
  lambda=lambda,
  decisions=decisions,
  lambda_MLE=CRRAMLEFit@fullcoef[1],
  lambda_MLEsd=CRRAMLEFit@vcov[1,1] %>% sqrt() %>% as.numeric(),
  lambda_B=sumCRRABayes["lambda","mean"],
  lambda_Bsd = sumCRRABayes["lambda","sd"],
  r_MLE=CRRAMLEFit@fullcoef[2],
  r_MLEsd=CRRAMLEFit@vcov[2,2] %>% sqrt() %>% as.numeric(),
  r_B=sumCRRABayes["r","mean"],
  r_Bsd = sumCRRABayes["r","sd"],
  Btime = CRRABayesFitTime,
  MLEtime = CRRAMLEFitTime
)
SimResults<-rbind(SimResults,
                  addThis)

addThis %>% print()
}
}

saveRDS(SimResults,file="ExampleStan/SimResults.rds")

SimResults<-readRDS("ExampleStan/SimResults.rds")

(
  ggplot(SimResults,aes(x=lambda_B,y=lambda_MLE))
  +theme_bw()
  +geom_point(alpha=0.2)
  +facet_grid(ifelse((r_B>0),"CRRA","Risk neutral")~decisions)
  +geom_abline(slope=1,intercept=0,linetype="dashed")
  +xlab("Bayesian estimate (posterior mean)")
  +ylab("Maximum likelihood estimate")
)
(
ggplot(SimResults,aes(x=lambda_Bsd,y=lambda_MLEsd))
+theme_bw()
+geom_point(alpha=0.2)
+facet_grid(ifelse((r_B>0),"CRRA","Risk neutral")~decisions)
+geom_abline(slope=1,intercept=0,linetype="dashed")
+xlab("Bayesian posterior standard deviation")
+ylab("Maximum likelihood standard error")
)


(
  ggplot(SimResults %>% filter(r_MLE>=0),aes(x=r_B,y=r_MLE))
  +geom_point(alpha=0.2)
  +theme_bw()
  +facet_wrap(~decisions)
  +geom_abline(slope=1,intercept=0,linetype="dashed")
  +xlab("Bayesian estimate (posterior mean)")
  +ylab("Maximum likelihood estimate")
)
