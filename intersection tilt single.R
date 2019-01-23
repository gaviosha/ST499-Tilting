# This function computes entry (v1,v2)
# of the tilted correlation matrix for X 
# using method 3.1 
# given threshold r

intersection_tilt_single <- function(v1, v2, X, r) {
  
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
  Z <- matrix(rep(1,nrow(X)), ncol = 1)
  for (i in 1:ncol(X)) {
    if (abs(t(X[,i]) %*% x1) > r 
        & abs(t(X[,i]) %*% x2) > r){
      Z <- cbind(Z, X[,i])
    }
  }
  
  #  if highly associated varaibles exist 
  #  find tilted correlation 
  if (ncol(Z) > 1){
    
    Z <- Z[,-1]
    P <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
    M <- diag(1,nrow(X)) - P
    
    #  transfrom x1 
    x1_str <- M %*% x1
    a <- t(x1) %*% M %*% x1 
    
    #  transfrom x2 
    x2_str <- M %*% x2
    b <- t(x2) %*% M %*% x2
    
    tilt <- (sqrt(a*b))^{-1} * (t(x1_str) %*% x2_str)
    
    # else return sraight correlation
  } else {
    tilt <- t(x1) %*% x2
  }
  
  return(tilt)
}