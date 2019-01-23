# This function computes entry (v1,v2)
# of the tilted correlation matrix for X 
# using method 3.3 
# given threshold r

average_tilt_single <- function(v1, v2, X, r) {
  
  #  transform data to have mean 0 and norm 1 
  X <- as.matrix(X)
  for (i in 1:ncol(X)) {
    X[,i] <- X[,i] - mean(X[,i])
    X[,i] <- X[,i]/sqrt(sum(X[,i]^{2}))
  }
  
  #  extract vectors and conditioning matrix
  x1 <- X[,v1]
  x2 <- X[,v2]
  X <- X[,-c(v1,v2)]
  
  # define empty conditioning set
  # find highly associated variables
  Z1 <- matrix(rep(1,nrow(X)), ncol = 1)
  Z2 <- matrix(rep(1,nrow(X)), ncol = 1)
  
  for (i in 1:ncol(X)) {
    if (abs(t(X[,i]) %*% x1) > r){
      Z1 <- cbind(Z1, X[,i])
    }
    if (abs(t(X[,i]) %*% x2) > r){
      Z2 <- cbind(Z2, X[,i])
    }
  }
  
  #  tilt x1  
  if (ncol(Z1) > 1){
    
    Z1 <- Z1[,-1]
    P1 <- Z1 %*% solve(t(Z1) %*% Z1) %*% t(Z1)
    M1 <- diag(1,nrow(X)) - P1
    
    x1_str <- M1 %*% x1
    a1 <- t(x1) %*% M1 %*% x1 
    
    tilt1 <- (t(x1_str) %*% x2)/a1
    
  } else {
    tilt1 <- t(x1) %*% x2
  }
  
  #  tilt x2  
  if (ncol(Z2) > 1){
    
    Z2 <- Z2[,-1]
    P2 <- Z2 %*% solve(t(Z2) %*% Z2) %*% t(Z2)
    M2 <- diag(1,nrow(X)) - P2
    
    x2_str <- M2 %*% x2
    a2 <- t(x2) %*% M2 %*% x2 
    
    tilt2 <- (t(x2_str) %*% x1)/a2
    
  } else {
    tilt2 <- t(x2) %*% x1
  }
  
  return(0.5*(tilt1 + tilt2))
}