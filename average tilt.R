# This function computes the complete tilted correlation matrix 
# using method 3.3 
# for a matrix X, using threshold r

average_tilt <- function(X, r) {
  
  #  transform data to have mean 0 and norm 1 
  X <- as.matrix(X)
  for (i in 1:ncol(X)) {
    X[,i] <- X[,i] - mean(X[,i])
    X[,i] <- X[,i]/sqrt(sum(X[,i]^{2}))
  }
  
  # produce tilted corelation matrix 
  tilt <- diag(1, ncol(X))
  
  for (i in 1:ncol(X)){
    for (j in 1:ncol(X)){
      if (i != j){
        
        # define new (x,y) variables 
        x <- X[,j]
        y <- X[,i]
        Xk <- X[,-c(i,j)]
        Z <- matrix(rep(1, nrow(X)), ncol = 1)
        
        # find highly associated variables 
        for (k in 1:ncol(Xk)) {
          if (abs(t(x) %*% Xk[,k]) > r) {
            Z <- cbind(Z, Xk[,k])
          }
        }
        
        if (ncol(Z) > 1) {
          
          #  if highly associated varaibles exist 
          #  define projection matrix 
          if (det(solve(t(Z) %*% Z)) == 0) stop("choose a higher threshold (r), 
                                                the projection matrix is singular!")
          Z <- Z[,-1]
          P <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
          M <- diag(1,nrow(X)) - P
          
          # project out highly associated
          x_str <- M %*% x
          a <- t(x) %*% M %*% x
          
          tilt[i,j] = (t(x_str) %*% y)/a
          
        } else {
          
          # If the matrix of associated regressors is empty... 
          # ...use straight correlation 
          tilt[i,j] = t(x) %*% y 
        } 
      }         
    }
  }
  
  # find average of correlation coeficients 
  for (i in 1:ncol(X)) {
    for (j in 1:ncol(X)) {
      if (i >= j) {
        tilt[i,j] <- 0.5*(tilt[i,j]+tilt[j,i])
        tilt[j,i] <- tilt[i,j]
      }
    }
  }
  return(tilt)
}