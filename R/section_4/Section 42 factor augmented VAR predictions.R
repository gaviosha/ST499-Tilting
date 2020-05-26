# load clean data and set X and y as desired 

dat <- read.csv("clean UK macro data.csv") 
dat <- dat[1:83,] # remove NaNs 
PCA <- princomp(dat[,2:27], cor = T) # extract principal components 
dat <- cbind(dat, PCA$scores[,1:3])
dat <- as.matrix(dat) # make numeric 
dat[,c(1,28:30)] <- dat[,c(1,28:30)]/100 # rescale CPI

Y <- cbind(dat[1:82,1],dat[2:83,c(12,26,109:111)])

# make variables to catch results 
pred.CPI <- NULL; MSE.CPI <- NULL; MAE.CPI <- NULL; size.CPI <- NULL
pred.XUQ <- NULL; MSE.XUQ <- NULL; MAE.XUQ <- NULL; size.XUQ <- NULL
pred.ABMI <- NULL; MSE.ABMI <- NULL; MAE.ABMI <- NULL; size.ABMI <- NULL

for (t in 1:32) {
  cat('just working on:',t,'\n')
  # set number of lags to 3 
  var.lags <- 3
  # estimate model and forecast 
  mod <- VAR(Y[t:(49+t),], p = var.lags, type = 'const')
  pred <- predict(mod, n.ahead = 1)$fcst
  pred <- as.numeric(as.data.frame(pred)[,c(1,5,9)])
  
  # fill in values for CPI
  pred.CPI <- c(pred.CPI,pred[1])
  MSE.CPI <- c(MSE.CPI,(pred[1]-Y[(t+50),1])^2)
  MAE.CPI <- c(MAE.CPI,abs(pred[1]-Y[(t+50),1]))
  size.CPI <- c(size.CPI, 3*var.lags)
  
  # fill in values for XUQABK67
  pred.XUQ <- c(pred.XUQ,pred[2])
  MSE.XUQ <- c(MSE.XUQ,(pred[2]-Y[(t+50),2])^2)
  MAE.XUQ <- c(MAE.XUQ,abs(pred[2]-Y[(t+50),2]))
  size.XUQ <- c(size.XUQ,3*var.lags)
  
  # fill in values for ABMI
  pred.ABMI <- c(pred.ABMI, pred[3])
  MSE.ABMI <- c(MSE.ABMI,(pred[3]-Y[(t+50),3])^2)
  MAE.ABMI <- c(MAE.ABMI,abs(pred[3]-Y[(t+50),3]))
  size.ABMI <- c(size.ABMI, 3*var.lags)
  
}

# put results in a data frame 
VAR.CPI.results <- data.frame(pred.CPI, MSE.CPI, MAE.CPI, size.CPI)
VAR.XUQ.results <- data.frame(pred.XUQ, MSE.XUQ, MAE.XUQ, size.XUQ)
VAR.ABMI.results <- data.frame(pred.ABMI, MSE.ABMI, MAE.ABMI, size.ABMI)

# summarize results
apply(VAR.CPI.results,2,mean)
apply(VAR.XUQ.results,2,mean)
apply(VAR.ABMI.results,2,mean)