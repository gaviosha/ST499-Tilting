#############################
# Tilting Simulations Setup #
#############################

# load required packages
require('tilting') # initial tilting package
require('tsDyn') # to simulate from stationary VAR
require('robustHD') # to standardize data
require('MASS'); require('mvtnorm') # to simulate multivariate normals and multivariate student-t
require('glmnet') # for lasso 
require('TAR') # for threshold autoregression
require('foreach'); require('doParallel') # for paralelisation
require('boot') # for stationary bootstrap 
require('TSA') # for simulating ARCH errors

##################################################################################################

# threshold function 
threshold <- function(X){
  
  # set base variables 
  n <- nrow(X)
  p <- ncol(X)
  FDR <- min(0.05, 1/sqrt(p))
  
  # get design matrix correlations
  C<-t(X)%*%X/n
  corr<-abs(C[upper.tri(C)])
  sc<-sort(corr, decreasing=FALSE)
  d<-p*(p-1)/2
  
  # make reference correlations 
  # generated from multivariate normal 
  D<-rmvnorm(n, sigma=diag(1, p))
  D<-t(t(D)/col.norm(D))
  C1<- t(D)%*%D
  ref<-abs(D[upper.tri(C1)])
  
  # find largest p-val satifsying test 
  i<-1
  while(i<=d){
    c<-sc[i]
    prob<-length(which(ref>c))/d
    if(prob<=(d-i+1)/d*FDR){
      return(c)
    }
    i<-i+1
  }
  
  # return associated correlation
  return(c)
}

##################################################################################################

# tilted corelation function 

tilt <- function(x,y,X, thr) {
  
  x <- standardize(x)
  y <- standardize(y)
  X <- apply(X,2, standardize)
  
  # find sample size and number of vars 
  n <- length(x)
  p <- ncol(X)
  
  # compute correlations between x and X 
  # keep columns of design matrix with correlations greater than pi 
  to_keep <- which(abs(cor(x,X)) > thr)
  to_keep <- to_keep[1:min(length(to_keep), sqrt(p))] # trim
  if (length(to_keep) == 0) {
    return(cor(x,y))
  }
  X <- X[,to_keep]
  
  # compute tilted correlation using re-scaling method 1 
  M = diag(1,n) - X%*%solve(t(X)%*%X)%*%t(X)
  z = M%*%x
  t_cor = (t(z)%*%y)/(t(x)%*%M%*%x)
  
  # return tilted correlation 
  return(t_cor)
}

##################################################################################################

# finsd tilted correlation with GLS-type weighting 

efficient.tilt <- function(x,y,X,thr) {
  
  # find sample size and number of vars 
  n <- length(x)
  p <- ncol(X)
  
  # compute correlations between x and X 
  # keep columns of design matrix with correlations greater than pi 
  to_keep <- which(abs(cor(x,X)) > thr)
  to_keep <- to_keep[1:min(length(to_keep), sqrt(p))] # trim
  
  # if no variables enter conditioning set return straigh correlation
  if (length(to_keep) == 0) {
    return(cor(x,y))
  }
  
  # else find efficient tilted correlation 
  res <- lm(x~X[,to_keep])$residuals # obtain residuals from regression on correlated predictors 
  AR <- ar(res) # fit AR(p) model to series of residuals 
  # check order of AR model
  if (AR$order == 0) {
    L <- diag(rep(1,n))*(1/sqrt(AR$var.pred))
  } else {
    S <- toeplitz(tacvfARMA(phi = AR$ar, maxLag = (n-1), sigma2 = AR$var.pred)) #for ACV-matrix 
    L <- solve(t(chol(S))) # inverse of colesky
  }
  
  # transform predictors and standardize 
  X <- L%*%X[,to_keep]; X <- apply(X,2, standardize)
  x <- L%*%x; x <- standardize(x)
  y <- L%*%y; y <- standardize(y)
  
  # compute tilted correlation using re-scaling method 1 
  M = diag(1,n) - X%*%solve(t(X)%*%X)%*%t(X)
  z = M%*%x
  t_cor = (t(z)%*%y)/(t(x)%*%M%*%x)
  
  # return tilted correlation 
  return(t_cor)
  
}

##################################################################################################

# finds tilted correlation with stable conditioning set 

stable.tilt <- function(x,y,X,thr) {
  
  # things we need 
  p <- ncol(X)
  n <- length(y)
  thr.seq <- seq(from = thr,to = 3*thr, length.out = 20)
  W <- matrix(rep(0,20*p),nrow = 20)
  
  for (i in 1:length(thr.seq)) {
    X.boot <- tsboot(tseries = data.frame(x,X), 
                     statistic = corr.extract, 
                     R = 50, 
                     sim = 'geom', 
                     l = 25)$t
    X.boot <- ifelse(X.boot>thr.seq[i],1,0) # which boostrap correlations exceed threshold 
    W[i,] <- apply(X.boot,2,mean) #  selection probability at threshold 
  }
  
  to_keep <- which(apply(W, 2, max) > 0.6)
  to_keep <- to_keep[1:min(length(to_keep), sqrt(p))] # trim
  
  # if no variables are selected by stability method return straigh correlation
  if (is.integer0(to_keep)) {
    return(cor(x,y))
  }
  
  # else compute tilted correlation 
  # transform predictors and standardize 
  X <- X[,to_keep]
  #X <- apply(X,2, standardize)
  x <- standardize(x)
  y <- standardize(y)
  
  # compute tilted correlation using re-scaling method 1 
  M = diag(1,n) - X%*%solve(t(X)%*%X)%*%t(X)
  z = M%*%x
  t_cor = (t(z)%*%y)/(t(x)%*%M%*%x)
  
  # return tilted correlation 
  return(t_cor)
  
}

##################################################################################################

# geometric series function 
GeometricSequence <- function(length, initial.value = 1, discount.factor = .5){
  stopifnot(is.numeric(length), length(length) == 1, length > 
              0, length == as.integer(length))
  stopifnot(is.numeric(initial.value), length(initial.value) == 
              1, initial.value != 0)
  stopifnot(is.numeric(discount.factor), length(discount.factor) == 
              1, discount.factor != 0)
  return(initial.value * discount.factor^(0:(length - 1)))
}

##################################################################################################

# theoretical autocovaraince of ARMA function 
# taken from repository for 'ltsa' package 
`tacvfARMA` <-
  function(phi = numeric(0), theta = numeric(0), maxLag = 1, sigma2=1)
  {
    stopifnot(maxLag >= 0, sigma2>0)
    #check for stationarity
    p <- length(phi)
    if(p>0){
      pi=numeric(p)
      phik <- phi
      for (k in 1:p){
        L <- p+1-k
        a <- phik[L]
        pi[L+1-k] <- a
        phikp1 <- phik[-L]
        if(is.na(a) || abs(a)==1) 
          stop("error: roots outside unit circle -- nonstationary/noncausal model")
        phik <- (phikp1+a*rev(phikp1))/(1-a^2)
      }
      if (!all(abs(pi)<1)) 
        stop("error: roots outside unit circle -- nonstationary/noncausal model")
    }
    #model is stationary, compute acvf
    q <- length(theta)
    if(max(p, q) == 0) 
      return(c(sigma2, numeric(maxLagp1))) #white noise case
    maxLagp1 <- maxLag + 1
    r <- max(p, q) + 1
    b <- numeric(r)
    C <- numeric(q + 1)
    C[1] <- 1
    theta2 <- c(-1, theta)
    phi2 <- numeric(3 * r)
    phi2[r] <- -1
    if(p > 0) 
      phi2[r + 1:p] <- phi
    if(q > 0) 
      for(k in 1:q) {
        C[k + 1] <-  - theta[k]
        if(p > 0) 
          for(i in 1:min(p, k)) 
            C[k + 1] <- C[k + 1] + phi[i] * C[k + 1 - i]
      }
    for(k in 0:q) 
      for(i in k:q) 
        b[k + 1] <- b[k + 1] - theta2[i + 1] * C[i - k + 1]
    if(p == 0) 
      g <- c(b, numeric(maxLagp1))[1:maxLagp1]
    else {
      a <- matrix(numeric(r^2), ncol = r)
      for(i in 1:r) 
        for(j in 1:r) {
          if(j == 1) 
            a[i, j] <- phi2[r + i - 1]
          else 
            a[i, j] <- phi2[r + i - j] + phi2[r + i + j - 2]
        }
      g <- solve(a,  - b)
      if(length(g) <= maxLag) {
        g <- c(g, numeric(maxLag - r))
        for(i in (r + 1):maxLagp1) 
          g[i] <- phi %*% g[i - 1:p]
      } 
    }
    sigma2*g[1:maxLagp1]
  }

##################################################################################################


# function for obtaining active set of given size using tilted correlations 
tilt.it.normal <- function(X,y,max.size, thr = threshold(X)){
  
  # standardize the data, extract dimensions 
  X <- apply(X,2,standardize)
  y <- as.matrix(standardize(y))
  n <- nrow(X)
  p <- ncol(X)
  
  # define objects for screening 
  Z <- X; X <- as.matrix(X) # design matrix 
  z <- y # residual 
  active <- NULL # active set
  I <- diag(rep(1,n))
  
  
  # itterate until derired size reached 
  while (length(active) < max.size){
    print(length(active))
    # first itteration: don't condition
    if (length(active)==0){
      k <- which.max(abs(cor(z,Z)))
      Ck <- which(abs(cor(X[,k],X)) > thr)
      
      # condition on next itterations  
    } else {
      k <- which.max(abs(cor(z,Z[,-active])))
      k <- setdiff(1:p,active)[k]
      Ck <- which(abs(cor(X[,k],X[,-active])) > thr)
      Ck <- setdiff(1:p,active)[Ck]
    }
    
    # If conditioning set only contains k 
    if (length(Ck)==1){
      # add k to active set
      active <- c(active,k)  
      # project out new active set 
      M <- as.matrix(X[,active])
      proj <- M%*%solve(t(M)%*%M)%*%t(M)
      z <- (I-proj)%*%y
      Z <- (I-proj)%*%as.matrix(X)
      Z <- as.data.frame(Z)
      # re-normalize columns of Z 
      Z <- apply(Z,2,standardize)
      
      # if conditioning set contains two ore more   
    } else {
      # find tilted correlation for each member of the set 
      t.cors <- c()
      for (j in Ck) {
        tc <- tilt(as.matrix(X[,j]),
                   z,
                   as.matrix(X[,-j]), 
                   thr)
        t.cors <- c(t.cors, tc)
      }
      # choose predictor with largest tilted 
      # correlation as k 
      k <- Ck[which.max(abs(t.cors))]
      # add k to active set
      active <- c(active,k)  
      # project out new active set
      M <- as.matrix(X[,active])
      # check if the matrix is singular first 
      if (abs(det(t(M)%*%M)) > 1.0e-3){
        proj <- M%*%solve(t(M)%*%M)%*%t(M)
        z <- (I-proj)%*%y
        Z <- (I-proj)%*%as.matrix(X)
      } else {
        e <- matrix(rnorm(ncol(M)*ncol(M), sd = 1.0e-3), 
                    nrow = ncol(M), 
                    ncol = ncol(M))
        proj <- M%*%solve(t(M)%*%M + e)%*%t(M)
        z <- (I-proj)%*%y
        Z <- (I-proj)%*%as.matrix(X)
      }
      Z <- as.data.frame(Z)
      # re-normalize columns of Z 
      Z <- apply(Z,2,standardize)
    } 
  }
  return(active)
}

##################################################################################################

# function for obtainig active set of a given size using efficient tilted correlations 
tilt.it.efficient <- function(X,y,max.size, thr = threshold(X)){
  
  # standardize the data, extract dimensions 
  X <- apply(X,2,standardize)
  y <- as.matrix(standardize(y))
  n <- nrow(X)
  p <- ncol(X)
  
  # define objects for screening 
  Z <- X; X <- as.matrix(X) # design matrix 
  z <- y # residual 
  active <- NULL # active set
  I <- diag(rep(1,n))
  
  
  # itterate until derired size reached 
  while (length(active) < max.size){
    print(length(active))
    # first itteration: don't condition
    if (length(active)==0){
      k <- which.max(abs(cor(z,Z)))
      Ck <- which(abs(cor(X[,k],X)) > thr)
      
      # condition on next itterations  
    } else {
      k <- which.max(abs(cor(z,Z[,-active])))
      k <- setdiff(1:p,active)[k]
      Ck <- which(abs(cor(X[,k],X[,-active])) > thr)
      Ck <- setdiff(1:p,active)[Ck]
    }
    
    # If conditioning set only contains k 
    if (length(Ck)==1){
      # add k to active set
      active <- c(active,k)  
      # project out new active set 
      M <- as.matrix(X[,active])
      proj <- M%*%solve(t(M)%*%M)%*%t(M)
      z <- (I-proj)%*%y
      Z <- (I-proj)%*%as.matrix(X)
      Z <- as.data.frame(Z)
      # re-normalize columns of Z 
      Z <- apply(Z,2,standardize)
      
      # if conditioning set contains two ore more   
    } else {
      # find tilted correlation for each member of the set 
      t.cors <- c()
      for (j in Ck) {
        tc <- efficient.tilt(as.matrix(X[,j]),
                             z,
                             as.matrix(X[,-j]), 
                             thr)
        t.cors <- c(t.cors, tc)
      }
      # choose predictor with largest tilted 
      # correlation as k 
      k <- Ck[which.max(abs(t.cors))]
      # add k to active set
      active <- c(active,k)  
      # project out new active set
      M <- as.matrix(X[,active])
      # check if the matrix is singular first 
      if (abs(det(t(M)%*%M)) > 1.0e-3){
        proj <- M%*%solve(t(M)%*%M)%*%t(M)
        z <- (I-proj)%*%y
        Z <- (I-proj)%*%as.matrix(X)
      } else {
        e <- matrix(rnorm(ncol(M)*ncol(M), sd = 1.0e-3), 
                    nrow = ncol(M), 
                    ncol = ncol(M))
        proj <- M%*%solve(t(M)%*%M + e)%*%t(M)
        z <- (I-proj)%*%y
        Z <- (I-proj)%*%as.matrix(X)
      }
      Z <- as.data.frame(Z)
      # re-normalize columns of Z 
      Z <- apply(Z,2,standardize)
    } 
  }
  return(active)
}

##################################################################################################

# function for obtaining active set of given size using stability tilting

# requires function to identify correlated predictors 
corr.extract <- function(X){
  return(abs(cor(X[,1],X[,-1])))
}

# requires parallel processing
registerDoParallel(cores=7) # register 7 cores as slaves and 1 as master 
getDoParWorkers() # check cores are registered


tilt.it.stable <- function(X,y,max.size,thr = threshold(X)){
  
  # standardize the data, extract dimensions 
  X <- apply(X,2,standardize)
  y <- as.matrix(standardize(y))
  n <- nrow(X)
  p <- ncol(X)
  
  # define objects for screening 
  Z <- X; X <- as.matrix(X) # design matrix 
  z <- y # residual 
  active <- NULL # active set
  I <- diag(rep(1,n))
  
  
  # itterate until derired Fstabsize reached 
  while (length(active) < max.size){
    print(length(active))
    # first itteration: don't condition
    if (length(active)==0){
      k <- which.max(abs(cor(z,Z)))
      Ck <- which(abs(cor(X[,k],X)) > thr)
      
      # condition on next itterations  
    } else {
      k <- which.max(abs(cor(z,Z[,-active])))
      k <- setdiff(1:p,active)[k]
      Ck <- which(abs(cor(X[,k],X[,-active])) > thr)
      Ck <- setdiff(1:p,active)[Ck]
    }
    
    # If conditioning set only contains k 
    if (length(Ck)==1){
      # add k to active set
      active <- c(active,k)  
      # project out new active set 
      M <- as.matrix(X[,active])
      proj <- M%*%solve(t(M)%*%M)%*%t(M)
      z <- (I-proj)%*%y
      Z <- (I-proj)%*%as.matrix(X)
      Z <- as.data.frame(Z)
      # re-normalize columns of Z 
      Z <- apply(Z,2,standardize)
      
      # if conditioning set contains two ore more   
    } else {
      
      # find tilted correlation for each member of the set 
      t.cors <- foreach (j = Ck, .combine = 'c') %dopar% {
        stable.tilt(as.matrix(X[,j]),
                    z,
                    as.matrix(X[,-j]), 
                    thr)
      }
      # choose predictor with largest tilted 
      # correlation as k 
      k <- Ck[which.max(abs(t.cors))]
      # add k to active set
      active <- c(active,k)  
      # project out new active set
      M <- as.matrix(X[,active])
      # check if the matrix is singular first 
      if (abs(det(t(M)%*%M)) > 1.0e-3){
        proj <- M%*%solve(t(M)%*%M)%*%t(M)
        z <- (I-proj)%*%y
        Z <- (I-proj)%*%as.matrix(X)
      } else {
        e <- matrix(rnorm(ncol(M)*ncol(M), sd = 1.0e-3), 
                    nrow = ncol(M), 
                    ncol = ncol(M))
        proj <- M%*%solve(t(M)%*%M + e)%*%t(M)
        z <- (I-proj)%*%y
        Z <- (I-proj)%*%as.matrix(X)
      }
      Z <- as.data.frame(Z)
      # re-normalize columns of Z 
      Z <- apply(Z,2,standardize)
    } 
  }
  return(active)
}


##################################################################################################

# function for obtaining active set of given size using stability tilting
# N.B.: screens tilted correlations of 100 most correlated predictors only! 

tilt.it.small.n.stable <- function(X,y,max.size,thr = threshold(X)){
  
  # standardize the data, extract dimensions 
  X <- apply(X,2,standardize)
  y <- as.matrix(standardize(y))
  n <- nrow(X)
  p <- ncol(X)
  
  # define objects for screening 
  Z <- X; X <- as.matrix(X) # design matrix 
  z <- y # residual 
  active <- NULL # active set
  I <- diag(rep(1,n))
  
  
  # itterate until derired Fstabsize reached 
  while (length(active) < max.size){
    print(length(active))
    # first itteration: don't condition
    if (length(active)==0){
      k <- which.max(abs(cor(z,Z)))
      Ck <- which(rank(abs(cor(X[,k],X))) >= 900)
      
      # condition on next itterations  
    } else {
      k <- which.max(abs(cor(z,Z[,-active])))
      k <- setdiff(1:p,active)[k]
      Ck <- which(rank(abs(cor(X[,k],X[,-active]))) >= 900)
      Ck <- setdiff(1:p,active)[Ck]
    }
    
    # If conditioning set only contains k 
    if (length(Ck)==1){
      # add k to active set
      active <- c(active,k)  
      # project out new active set 
      M <- as.matrix(X[,active])
      proj <- M%*%solve(t(M)%*%M)%*%t(M)
      z <- (I-proj)%*%y
      Z <- (I-proj)%*%as.matrix(X)
      Z <- as.data.frame(Z)
      # re-normalize columns of Z 
      Z <- apply(Z,2,standardize)
      
      # if conditioning set contains two ore more   
    } else {
      
      # find tilted correlation for each member of the set 
      t.cors <- foreach (j = Ck, .combine = 'c') %dopar% {
        stable.tilt(as.matrix(X[,j]),
                    z,
                    as.matrix(X[,-j]), 
                    thr)
      }
      # choose predictor with largest tilted 
      # correlation as k 
      k <- Ck[which.max(abs(t.cors))]
      # add k to active set
      active <- c(active,k)  
      # project out new active set
      M <- as.matrix(X[,active])
      # check if the matrix is singular first 
      if (abs(det(t(M)%*%M)) > 1.0e-3){
        proj <- M%*%solve(t(M)%*%M)%*%t(M)
        z <- (I-proj)%*%y
        Z <- (I-proj)%*%as.matrix(X)
      } else {
        e <- matrix(rnorm(ncol(M)*ncol(M), sd = 1.0e-3), 
                    nrow = ncol(M), 
                    ncol = ncol(M))
        proj <- M%*%solve(t(M)%*%M + e)%*%t(M)
        z <- (I-proj)%*%y
        Z <- (I-proj)%*%as.matrix(X)
      }
      Z <- as.data.frame(Z)
      # re-normalize columns of Z 
      Z <- apply(Z,2,standardize)
    } 
  }
  return(active)
}

##################################################################################################

# functions for ROC analysis 

make.TPR <- function(mod.seq){
  required <- c(200,400,600,800,1000)
  TPR <- 0 
  for (i in 1:length(mod.seq)){
    tp <- length(which(mod.seq[1:i] %in% required))/5
    TPR <- c(TPR, tp)
  }
  return(TPR)
}

make.FPR <- function(mod.seq){
  required <- c(200,400,600,800,1000)
  FPR <- 0 
  for (i in 1:length(mod.seq)){
    fp <- length(mod.seq[1:i][-which(mod.seq[1:i] %in% required)])/995
    FPR <- c(FPR, fp)
  }
  return(FPR)
}

is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

