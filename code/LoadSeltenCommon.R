# LoadSeltenCommon.R
# Common starting point for all scripts so I'm working with the same cleaned data
#install.packages("tidyr")
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
options(mc.cores = 4) # Raise if you have lots of RAM. Lower if you don't
library(ggthemes)
library(ggplot2)
library(latex2exp)
library(hexbin)
library(tidyverse)

dataIndividual<-(read_delim("2x2_all_data.csv",delim=";")
       %>% mutate(role = ifelse(!is.na(Left),"column","row"))
       %>% filter(round>=101) # focus on last half of decisions
       %>% select(game,observation,round,player,role,LeftUp)
       %>% group_by(game,observation,player,role)
       %>% summarize(LeftUp = sum(LeftUp),count=n())
       %>% pivot_wider(values_from = LeftUp,names_from=role)
       
)

Up_i<-dataIndividual %>% ungroup()%>% filter(!is.na(row)) %>% select(game,row,count)
Left_i<-dataIndividual %>% ungroup() %>% filter(!is.na(column)) %>% select(game,column,count)


data<-(read_delim("2x2_all_data.csv",delim=";")
       %>% mutate(role = ifelse(!is.na(Left),"column","row"))
       %>% filter(round>=101) # focus on last half of decisions
  %>% select(game,observation,round,player,role,LeftUp)  
  %>% group_by(game,observation,role)
  %>% summarize(LeftUp = sum(LeftUp),count = n())
  %>% pivot_wider(values_from = LeftUp,names_from=role)
)


#data <- read_excel("../data/QREMeta.xlsx",sheet="Data",col_names = TRUE)
#data <- read_excel("data/Selten2008.xlsx",sheet="Data",col_names = TRUE)
# set ==1 to use the data, otherwise I just check the prior specifications
#run_estimation<-1

#-------------------------------------------------------------------
# get the data into a more JB-friendly format
#-------------------------------------------------------------------



UDLRcount<-cbind(data$row,data$count-data$row,data$column,data$count-data$column)

A <- rbind(c(1,-1,0,0),c(0,0,1,-1))

np=30


lxgrid<-seq(0,0.5,length=np)/(1-seq(0,0.5,length=np))

SC2008payoffs<-list()
SC2008payoffs[[1]]<-list(
  Urow = c(10,0,9,10) %>% matrix(nrow=2) %>% t(),
  Ucol = c(8,18,9,8) %>% matrix(nrow=2) 
)
SC2008payoffs[[2]]<-list(
  Urow = c(9,0,6,8) %>% matrix(nrow=2) %>% t(),
  Ucol = c(4,13,7,5) %>% matrix(nrow=2) 
)
SC2008payoffs[[3]]<-list(
  Urow = c(8,0,7,10) %>% matrix(nrow=2)%>% t(),
  Ucol = c(6,14,7,4) %>% matrix(nrow=2) 
)
SC2008payoffs[[4]]<-list(
  Urow = c(7,0,5,9) %>% matrix(nrow=2) %>% t(),
  Ucol = c(4,11,6,2) %>% matrix(nrow=2) 
)
SC2008payoffs[[5]]<-list(
  Urow = c(7,0,4,8) %>% matrix(nrow=2)%>% t(),
  Ucol = c(2,9,5,1) %>% matrix(nrow=2) 
)
SC2008payoffs[[6]]<-list(
  Urow = c(7,1,3,8) %>% matrix(nrow=2) %>% t(),
  Ucol = c(1,7,5,0) %>% matrix(nrow=2)
)
SC2008payoffs[[7]]<-list(
  Urow = c(10,4,9,14) %>% matrix(nrow=2)%>% t(),
  Ucol = c(12,22,9,8) %>% matrix(nrow=2) 
)
SC2008payoffs[[8]]<-list(
  Urow = c(9,3,6,11) %>% matrix(nrow=2) %>% t(),
  Ucol = c(7,16,7,5) %>% matrix(nrow=2)
)
SC2008payoffs[[9]]<-list(
  Urow = c(8,3,7,13) %>% matrix(nrow=2) %>% t(),
  Ucol = c(9,17,7,4) %>% matrix(nrow=2)
)
SC2008payoffs[[10]]<-list(
  Urow = c(7,2,5,11) %>% matrix(nrow=2)%>% t(),
  Ucol = c(6,13,6,2) %>% matrix(nrow=2) 
)
SC2008payoffs[[11]]<-list(
  Urow = c(7,2,4,10) %>% matrix(nrow=2)%>% t(),
  Ucol = c(4,11,5,1) %>% matrix(nrow=2) 
)
SC2008payoffs[[12]]<-list(
  Urow = c(7,3,3,10) %>% matrix(nrow=2)%>% t(),
  Ucol = c(3,9,5,0) %>% matrix(nrow=2) 
)

Payoffs<-c()
for (gg in SC2008payoffs) {
  Payoffs<-rbind(Payoffs,
    c(as.vector(gg[[1]] %>% t()),as.vector(gg[[2]] %>% t()))
      )
}

d<-list(
  Payoffs = Payoffs,
  nc=20,
  ftol=1e-5,
  A=A,
  
  UDLRcount=UDLRcount,
  GameID=data$game,
  ncohorts=dim(UDLRcount)[1],
  ngames=length(unique(data$game)),
  
  prior_lambda=c(-1.52,1.41),
  prior_r=c(log(0.5),0.5), 
  
  prior_taulambdaGame = c(log(0.3),0.2),
  prior_taurGame = c(log(0.3),0.2),
  prior_taulambdaCohort = c(log(0.3),0.2),
  prior_taurCohort = c(log(0.3),0.2),
  prior_taulambdaBoth = c(0.5*log(0.3),0.2/sqrt(2)),
  prior_taurBoth = c(0.5*log(0.3),0.2/sqrt(2)),
  
  prior_lambda_sd_game = c(1,12),
  prior_lambda_sd_cohort = c(1,12),
  prior_lambda_sd_both = c(0.5,12),
  prior_r_sd_game = c(0.5,10),
  prior_r_sd_cohort = c(0.5,10),
  prior_r_sd_both = c(0.25,10),
  priorLogitN = c(log(0.1),0.1,3),
  prior_LNsd = 0.25, # Think about this one
  prior_lkj = 2, 
  prior_lnSD = c(0.1,3), # Gamma prior
  prior_kappa = 0.0047,
  
  #pNash=cbind(pNashU,pNashL),
  
  
  UseData=rep(1,length(data$game)),
  
  ZR = qnorm(seq(0.01,0.99,length=30)),
  nZ = 30,
  Up_i = Up_i$row,
  Left_i = Left_i$column,
  count_i = Up_i$count,
  game_i = Up_i$game,
  nplayers = dim(Up_i)[1],
  prior_HQREsd = 8.6
)

saveRDS(d,"SeltenChmura2008data.rds")


