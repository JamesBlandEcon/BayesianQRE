rm(list=ls()) 
library("rstan") # observe startup messages
rstan_options(auto_write = TRUE)
library(stringr)
library(tidyr)
library(dplyr)


recompile<-FALSE

recompileList<-c()

ModelFiles<-list.files(path="StanFiles",pattern="_Ind.stan")

ModelFiles<-ModelFiles[!grepl(".stani",ModelFiles)]

ModelList<-data.frame()
for (mm in ModelFiles) {
  m<-unlist(str_split(mm,".stan"))[1]
  
  tmp<-data.frame(t(unlist(str_split(m,"_"))))
  colnames(tmp)<-c("Solution concept","Risk","Heterogeneity","Correlation")
  tmp$StanFile <-mm
  
  ModelList<-rbind(ModelList,tmp)
}

ModelList<-ModelList %>% arrange(`Solution concept`,Heterogeneity,Risk)

print(ModelList)


if (!file.exists("Models.rds") | recompile) {
  MODELS<-list()
  saveRDS(MODELS,"Models.rds")
}

MODELS<-readRDS("Models.rds")



for (mm in ModelFiles) {
  if (length(MODELS)==0) {
    SM<-0
  }
  else {
    SM<-sum(names(MODELS)==mm)
  }
  if (SM==0 | any(mm==recompileList) ) {
  print(paste("Starting to compile model ",mm,"at",Sys.time()))
  MODELS[[mm]]=stan_model(paste0("StanFiles/",mm),
                          model_name=mm,
                          auto_write=TRUE,
                          allow_undefined=TRUE,
                          verbose=FALSE,
                          includes=c("\n"))
  saveRDS(MODELS,"Models.rds")
  print(paste("model ",mm,"compiled at",Sys.time()))
  
  }
  else {
    paste("no need to compile model ",mm,"again :)")
  }
}

source("EstimateAll.R")


