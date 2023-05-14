source("LoadSeltenCommon.R")

reestimate = FALSE


ReEstimateList<-c()

controlList<-list(adapt_delta=0.999,max_treedepth=20)

Models<-readRDS("Models.rds")

if (!file.exists("Estimates.rds") | reestimate) {
  ESTIMATES<-list()
  saveRDS(ESTIMATES,"Estimates.rds")
} else {
  ESTIMATES<-readRDS("Estimates.rds")
}


for (mm in names(Models)) {
  if ((ESTIMATES[[mm]] %>% is.null()) | any(mm==ReEstimateList)) {
    paste("Estimating model ",mm) %>% print()
    data<-d
    if (grepl("NBE",mm,fixed=TRUE)) {
      # Override prior lambda for NBE model
      data$prior_lambda<-c(-0.2221069, 0.9838724)
    } else if (grepl("Luce",mm,fixed=TRUE)) {
      data$prior_lambda<-c( 1.497461, 1.781139)
    }
    
    Fit<-sampling(
      Models[[mm]],
      data  = data,
      control = controlList,
      seed = 1234
    )
    
    ESTIMATES[[mm]]<-Fit
    saveRDS(ESTIMATES,"Estimates.rds")
  } else {
    paste("Skipping model ",mm) %>% print()
  }
}