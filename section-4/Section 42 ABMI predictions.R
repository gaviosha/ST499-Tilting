# load clean data and set X and y as desired 

dat <- read.csv("clean UK macro data.csv") 
dat <- dat[1:83,] # remove NaNs 
PCA <- princomp(dat[,2:27], cor = T) # extract principal components 
dat <- cbind(dat, PCA$scores[,1:3])
dat <- as.matrix(dat) # make numeric 
dat[,c(1,28:30)] <- dat[,c(1,28:30)]/100 # rescale CPI 
X <- dat[1:82,-c(1,26)] # assign design matrix 
y <- dat[2:83,26] # assign target

## calculations for AR(p) model 
##
# make variables to catch results 
pred <- NULL; MSE <- NULL; MAE <- NULL; size <- NULL

for (t in 1:32){
  # fit AR model with AIC 
  mod <- ar(y[t:(49+t)])
  # obtain prediction 
  fitted <- as.numeric(predict(mod, n.ahead = 1)$pred)
  # fill in metrics 
  pred <- c(pred, fitted)
  MSE <- c(MSE, (y[50+t]-fitted)^2)
  MAE <- c(MAE,abs(y[50+t]-fitted))
  size <- c(size, mod$order)
}

# view results
ar.results <- data.frame(y[51:82],pred,MSE,MAE,size)
View(ar.results)
apply(ar.results,2,mean)


## calculations for factor augmented AR model 
##
# make variables to catch results 
pred <- NULL; MSE <- NULL; MAE <- NULL; size <- NULL

for (t in 1:32){
  # select three lags and three factors  
  to_choose <- c(101:103,107:109) # choose first three lags and three factors
  cat('working on:',t,'for factro augmented AR model','\n')
  mod <- glmnet(X[t:(49+t),to_choose],y[t:(49+t)], alpha = 1, lambda = 0)
  # obtain prediction 
  fitted <- predict.glmnet(object = mod, newx = t(X[50+t,to_choose]), s = 0)
  # fill in metrics 
  pred <- c(pred, fitted)
  MSE <- c(MSE, (y[50+t]-fitted)^2)
  MAE <- c(MAE,abs(y[50+t]-fitted))
  size <- c(size, length(which(mod$beta != 0)))
}

# View results 
factor.ar.results <- data.frame(y[51:82],pred,MSE,MAE,size)
View(factor.ar.results)
apply(factor.ar.results,2,mean)

## calculations for normal tilting 
##
# make variables to catch results 
pred <- NULL; MSE <- NULL; MAE <- NULL; size <- NULL

for (t in 1:32){
  # select variables via tilting, fit with lasso 
  to_choose <- tilt.it.normal(X[t:(49+t),],y[t:(49+t)], max.size = 25)
  to_choose <- union(to_choose,27:29)
  cat('working on:',t,'for normal tilting','\n')
  opt.lambda <- cv.glmnet(X[t:(49+t),],y[t:(49+t)])
  mod <- glmnet(X[t:(49+t),to_choose],y[t:(49+t)], alpha = 1, lambda = opt.lambda$lambda.min)
  # obtain prediction 
  fitted <- predict.glmnet(object = mod, newx = t(X[50+t,to_choose]), s = opt.lambda$lambda.min)
  # fill in metrics 
  pred <- c(pred, fitted)
  MSE <- c(MSE, (y[50+t]-fitted)^2)
  MAE <- c(MAE,abs(y[50+t]-fitted))
  size <- c(size, length(which(mod$beta != 0)))
}

# View results
tilting.normal.results <- data.frame(y[51:82],pred,MSE,MAE,size)
View(tilting.normal.results)
apply(tilting.normal.results,2,mean)

## calculations for efficient tilting 
##
# make variables to catch results 
pred <- NULL; MSE <- NULL; MAE <- NULL; size <- NULL

for (t in 1:32){
  # select variables via tilting, fit with lasso 
  to_choose <- tilt.it.efficient(X[t:(49+t),],y[t:(49+t)], max.size = 25)
  cat('working on:',t,'for efficient tilting','\n')
  opt.lambda <- cv.glmnet(X[t:(49+t),],y[t:(49+t)])
  mod <- glmnet(X[t:(49+t),to_choose],y[t:(49+t)], alpha = 1, lambda = opt.lambda$lambda.min)
  # obtain prediction 
  fitted <- predict.glmnet(object = mod, newx = t(X[50+t,to_choose]), s = opt.lambda$lambda.min)
  # fill in metrics 
  pred <- c(pred, fitted)
  MSE <- c(MSE, (y[50+t]-fitted)^2)
  MAE <- c(MAE,abs(y[50+t]-fitted))
  size <- c(size, length(which(mod$beta != 0)))
}

# View results 
tilting.efficient.results <- data.frame(y[51:82],pred,MSE,MAE,size)
View(tilting.efficient.results)
apply(tilting.efficient.results,2,mean)

## calculations for stability tilting 
##
# make variables to catch results 
pred <- NULL; MSE <- NULL; MAE <- NULL; size <- NULL

for (t in 1:32){
  # select variables via tilting, fit with lasso 
  to_choose <- tilt.it.stable(X[t:(49+t),],y[t:(49+t)], max.size = 25)
  cat('working on:',t,'for efficient tilting','\n')
  opt.lambda <- cv.glmnet(X[t:(49+t),],y[t:(49+t)])
  mod <- glmnet(X[t:(49+t),to_choose],y[t:(49+t)], alpha = 1, lambda = opt.lambda$lambda.min)
  # obtain prediction 
  fitted <- predict.glmnet(object = mod, newx = t(X[50+t,to_choose]), s = opt.lambda$lambda.min)
  # fill in metrics 
  pred <- c(pred, fitted)
  MSE <- c(MSE, (y[50+t]-fitted)^2)
  MAE <- c(MAE,abs(y[50+t]-fitted))
  size <- c(size, length(which(mod$beta != 0)))
}

# View results
tilting.stable.results <- data.frame(y[51:82],pred,MSE,MAE,size)
View(tilting.stable.results)
apply(tilting.stable.results,2,mean)

# View summary of results 
apply(ar.results,2,mean)
apply(factor.ar.results,2,mean)
apply(tilting.normal.results,2,mean)
apply(tilting.efficient.results,2,mean)
apply(tilting.stable.results,2,mean)
