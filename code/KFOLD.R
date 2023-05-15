source("LoadSeltenCommon.R")
library(loo)

reestimate = FALSE


ReEstimateList<-c()

controlList<-list(adapt_delta=0.999,max_treedepth=20)

Models<-readRDS("Models.rds")

set.seed(42)

#folds<-ceiling((1:d$ncohorts)/18)[order(runif(d$ncohorts))]
folds<-d$GameID

if (!file.exists("kfold.rds") | reestimate) {
  LPLIST<-list()
  saveRDS(LPLIST,"kfold.rds")
} else {
  LPLIST<-readRDS("kfold.rds")
}



for (ff in unique(folds)) {
  dfold<-d
  dfold$UseData<-folds!=ff
  for (mm in names(Models)) {
    if (is.null(LPLIST[[mm]])|any(is.na(LPLIST[[mm]][,folds==ff]))) {
    paste("Estimating fold",ff, "model ",mm) %>% print()
      data<-dfold
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
      seed = 12345
    )
    lp<-extract_log_lik(Fit)
    if (length(LPLIST)==0 | is.null(LPLIST[[mm]])) {
      lpBlank<-NA*lp
      
        LPLIST[[mm]]<-lpBlank
      }
    LPLIST[[mm]][,folds==ff]<-lp[,folds==ff]
    saveRDS(LPLIST,file="kfold.rds")
    } else {
      paste("skipping fold",ff, "model ",mm,"- already done") %>% print()
    }
  }
}