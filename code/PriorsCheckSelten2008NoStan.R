# Calibrate the priors. This script does not use Stan, so will run faster 
# when there is no likelihood, but should produce the same results
rm(list=ls()) 
library("rstan") # observe startup messages
rstan_options(auto_write = TRUE)
library(foreign)
library(png)
library(grid)
library(xtable)
library("readxl")
#library(bridgesampling)
library(rootSolve)
options(mc.cores = min(c(parallel::detectCores(),4)))
library(ggthemes)
library(ggplot2)
library(latex2exp)
library(hexbin)
library(gganimate)
library(dplyr)
library(tidyr)
#library(directlabels)
#library(tweenr)
data <- read_excel("../data/Selten2008.xlsx",sheet="Data",col_names = TRUE)

GameData<-data %>% group_by(Game) %>% filter(Session==min(Session)) %>% select(-U,-L,-N,-Session)

GameList<-list()
for (gg in 1:dim(GameData)[1]) {
  ur<-GameData[gg,2:5]
  uc<-GameData[gg,6:9]
  
  Urow<-rbind(c(ur$UL_R,ur$UR_R),c(ur$DL_R,ur$DR_R))
  Ucol<-t(rbind(c(uc$UL_C,uc$UR_C),c(uc$DL_C,uc$DR_C)))
  
  GameList[[gg]]<-list(Urow,Ucol)
  
}

## Determine what r ~ lognormal(m,s) means at Nash
r<-seq(0.1,1.5,length=15)
CRRA<-data.frame()
for (gg in 1:length(GameList)) {
  Urow<-GameList[[gg]][[1]]
  Ucol<-GameList[[gg]][[1]]
  for (rr in r) {
    pL<-solve(Urow^rr,c(1,1))
    pL<-pL[1]/sum(pL)
    pU<-solve(Ucol^rr,c(1,1))
    pU<-pU[1]/sum(pU)
    CRRA<-data.frame(pU,pL) %>% mutate(r=rr,Game=gg) %>% rbind(CRRA) 
  }
}

plt<-( 
  ggplot(data=CRRA ,aes(x=pU,y=pL,color=r))
  +geom_path()
  +geom_text(data=CRRA %>% filter(abs(r-1)<0.01),aes(label="rn"))
  +facet_wrap(~Game)
  )
print(plt)


MakePlot<-0
L_ii = seq(0.001,0.999,length=1000)
LAMBDA    <- L_ii/(1-L_ii)


GameList <-unique(data$Game)
PayoffList<-matrix(0,length(GameList),8)

QRE<-c()

QRELabels<-c("lambda")

LX<-as.matrix(LAMBDA)
LXCols<-c("lambda")

for (gg in 1:max(GameList)) {
  print(gg)
  qre <-matrix(-1,length(LAMBDA),2)
  
  Ugg<-colMeans(data[data$Game==gg,6:13])
  Urow<-rbind(c(Ugg[1],Ugg[2]),c(Ugg[3],Ugg[4]))
  Ucol<-t(rbind(c(Ugg[5],Ugg[6]),c(Ugg[7],Ugg[8])))
  
  for (ll in 1:length(LAMBDA)) {
    l<-LAMBDA[ll]
    obj<-function(x) {
      pU <- 1.0/(1+exp(-x[1]))
      pL <- 1.0/(1+exp(-x[2]))
      F1 <- x[1]-l*((Urow[1,1]-Urow[2,1])*pL+(Urow[1,2]-Urow[2,2])*(1.0-pL));
      F2 <- x[2]-l*((Ucol[1,1]-Ucol[2,1])*pU+(Ucol[1,2]-Ucol[2,2])*(1.0-pU));
      c(F1=F1,F2=F2)
    }
    sol <- multiroot(f=obj,start(c(0,0)))
    qre[ll,1] <- 1./(1+exp(-sol$root[1]))
    qre[ll,2] <- 1./(1+exp(-sol$root[2]))
  }
  
  dP<-rowSums((qre[2:length(LAMBDA),]-qre[1:(length(LAMBDA)-1),])^2)
  x<-cumsum(dP)/sum(dP)
  x<-c(0,x)
  
  QRE<-rbind(QRE,cbind(matrix(gg,length(LAMBDA),1),LAMBDA,qre,x))
  
  LX<-cbind(LX,x)
  LXCols<-c(LXCols,paste("x",gg,sep=""))
  
}

## QRE plot
colnames(QRE)<-c("Game","Lambda","pU","pL","x")
QREdata<-data.frame(QRE)
saveRDS(QREdata,file="QRETrace.rds")
plt<-ggplot(QREdata,aes(x=pU,y=pL,group=Game))+geom_path()+theme(legend.position = "none")
plt<-plt+xlab("Pr[U]")+ylab("Pr[L]")+theme_bw()
ggsave(paste("../outputs/priorplots/SeltenChmura_QREtrace.png", sep=""), width = 5, height = 5)
plt

colnames(QRE)<-c("Game","lambda","pU","pL","x")
colnames(LX)<-LXCols

med<-c()
for (gg in 1:length(GameList)) {
  med<-rbind(med,which.min(abs(LX[,gg]-0.5)))
}
m<-log(median(LAMBDA[med]))

set.seed(42)
Z <- rnorm(10000)


## Objective function to minimize KS-statistic
KS<- function(xx) {
  mu <- xx[1]
  sigma <-xx[2]
  Lsim<-exp(mu+sigma*Z)
  ks<-c()
  for (gg in 1:length(GameList)) {
    finterp<-approxfun(LX[,1],LX[,gg+1],method="linear",yleft=0,yright=1)
    xsim<-finterp(Lsim)
    P<-ecdf(xsim)
    ks<-c(ks,max(abs(xsim-P(xsim))))
    
    
  }
  
  c(max(ks))
}


opt<-optim(c(m,1),KS,gr = NULL)







Xsim<-c()
ggTemp<-c()
Mtemp<-c()
Stemp<-c()


S<- round(opt$par[2]*c(0.2,0.6,1.0,1.8,2.5),2)
M<- round(opt$par[1]*c(0.2,0.6,1.0,1.4,1.8),2)

#S<- round(opt$par[2]*c(0.6,1.0,1.4),2)
#M<- round(opt$par[1]*c(0.6,1.0,1.4),2)



for (ss in 1:length(S)) {
for (mm in 1:length(M)) {
  s<-S[ss]
  m<-M[mm]
  Lsim<-exp(m+s*Z)
  for (gg in 1:max(GameList)) {
    finterp<-approxfun(LX[,1],LX[,gg+1],method="linear",yleft=0,yright=1)
    xx<-finterp(Lsim)
    Xsim<-rbind(Xsim,as.matrix(xx))
    ggTemp<-rbind(ggTemp,matrix(paste(gg,"",sep=""),length(Lsim),1))
    Mtemp<-rbind(Mtemp,matrix(m,length(Lsim),1))
    Stemp<-rbind(Stemp,matrix(s,length(Lsim),1))
    
    
    
  }
  print(c(ss,mm))
}
  
} 

df<-data.frame(Xsim)
df$M<-Mtemp
df$S<-Stemp
df$game<-ggTemp

if (MakePlot==1) {

plt<-ggplot(df,aes(Xsim,color=game, parse=TRUE))+geom_abline(slope=1, intercept=0,linetype=2)+theme_bw() +stat_ecdf(geom="step") + facet_grid(M ~ S, labeller = label_both)
plt<-plt+xlab(TeX("$x_g(\\lambda)$: Normalized distance from $\\lambda = 0$"))+ylab("Cumulative density")+theme(legend.position = "none")
ggsave(paste("../outputs/priorplots/SeltenChmura_Homogeleous_PriorCalibration.png", sep=""), width = 10, height = 10)
plt
}

## Determine Gamma prior on Sigma
z<-qnorm(0.75)
DX<-c()
XX<-c()
GG<-c()
SS<-c()
SIGMA<-c(0.25,0.5,0.75)
for (ss in 1:length(SIGMA))  {
  sig<-qgamma(SIGMA[ss], 1, 12, lower.tail = TRUE,log.p = FALSE)

  for (gg in 1:length(GameList)) {
    finterp<-approxfun(LX[,1],LX[,gg+1],method="linear",yleft=0,yright=1)
    dx<-finterp(LX[,1]+sig*z)-finterp(LX[,1]-sig*z)
    DX<-c(DX,dx)
    XX<-c(XX,as.vector(LX[,gg+1]))
    GG<-c(GG,as.vector(matrix(gg+SIGMA[ss],length(dx),1)))
    SS<-c(SS,as.vector(matrix(paste(SIGMA[ss]*100,"th",sep=""),length(dx),1)))
  }
  }


df<-data.frame(GG)
df$DX<-DX
df$XX<-XX
df$percentile<-SS

plt<-ggplot(df,aes(x=XX,y=DX,group=GG,color=percentile))+geom_line()+theme_bw()
plt<-plt+xlab(TeX("$x_g(\\bar{\\lambda})$: Normalized distance from $\\lambda = 0$"))+ylab("Interquartile Range Coverage")
plt
ggsave(paste("../outputs/priorplots/SeltenChmura_Hierarchical_PriorCalibration.png", sep=""), width = 5, height = 5)
stop()
# Equilibrium in expectation models

prior_K<-c(1,2) # initial guess

obj<-function(xx) {
  K<-exp(xx[1]+xx[2]*Z)
  Ksd<-1/sqrt(1+K)
  P<-ecdf(Ksd)
  max(abs(Ksd-P(Ksd)))
}

opt<-optim(prior_K,obj,gr = NULL)
prior_K<-opt$par

kSDsim<-1/sqrt(1+exp(prior_K[1]+prior_K[2]*Z))

df<-data.frame(kSDsim)
plt<-ggplot(data=df,aes(x=kSDsim))+geom_abline(slope=1, intercept=0,linetype=2)+theme_bw() +stat_ecdf(geom="step")
plt<-plt+xlab(TeX("$\\frac{1}{\\sqrt{1+\\kappa}}$"))+ylab("Cumulative density")
plt
ggsave(paste("../outputs/priorplots/SeltenChmura_K_PriorCalibration.png", sep=""), width = 5, height = 5)


d<-QRE[QRE[,2]== max(QRE[,2]),]

pNash<-unique(as.vector(d[,3:4]))

qre<-QRE[LAMBDA==LAMBDA[1] | LAMBDA==LAMBDA[length(LAMBDA)/2] | LAMBDA==LAMBDA[length(LAMBDA)],]

# Minimum K to ensure distribution of p is single-peaked
minK<-1/min(c(pNash,1-pNash))

x<-seq(0.01,0.99,0.01)
y<-dbeta(x,0.2 * minK,(1-0.2) * minK)
d<-data.frame(x,y)
d$k<-"minimum"
d$p<-"0.2"

tmp<-d
tmp$y<-dbeta(x,0.5 * minK,0.5 * minK)
tmp$p<-"0.5"
d<-rbind(d,tmp)

tmp$y<-dbeta(x,0.2 * 2*minK,(1-0.2) * 2*minK)
tmp$p<-"0.2"
tmp$k<-"2x minimum"
d<-rbind(d,tmp)

tmp<-d
tmp$y<-dbeta(x,0.5 * 2*minK,0.5 * 2*minK)
tmp$p<-"0.5"
tmp$k<-"2x minimum"
d<-rbind(d,tmp)

plt<-(
  ggplot(data=d,aes(x=x,y=y,color=p,linetype=k))
  + geom_line(size=1)
  +theme_bw()
  
)
plt

U<-runif(10000)
f<-function(a) {
sqrt(1/4/(minK/(1-U)^a+1))
}

target<-function(a){
  sdSim<-f(a)
  P<-ecdf(sdSim)(sdSim)
  max(abs(P-sdSim/sqrt(1/4/(minK+1))))
}

a<-optim(3,target,gr = NULL,method="Brent",lower=1,upper=100)


AA<-c(1,a$par,3)
d<-data.frame()
for (aa in 1:length(AA)) {
  sdSim<-f(AA[aa])
  P<-ecdf(sdSim)(sdSim)
  tmp<-data.frame(sdSim,P)
  tmp$a<-paste(round(AA[aa],1))
  d<-rbind(d,tmp)
  
}

plt<-(
  ggplot(data=d,aes(x=sdSim,y=P,linetype=a))+geom_line(size=1)
  +geom_abline(slope=1/sqrt(1/4/(minK+1)),intercept=0,color="red")
  +theme_bw()
  +scale_linetype(name="\U03B1")
  + theme(legend.position = c(0.8, 0.2))
  +xlab("Standard deviation of Beta distribution")+ylab("Prior cumulative density function")
)
plt
ggsave(paste("../outputs/priorplots/SeltenChmura_KBETA_PriorCalibration.png", sep=""), width = 5, height = 5)
