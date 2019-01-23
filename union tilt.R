# This function computes the complete tilted correlation matrix 
# using method 3.2 
# for a matrix X, using threshold r

union_tilt <- function(X,r) {
  
  #  transform data to have mean 0 and norm 1 
  X <- as.matrix(X)
  for (i in 1:ncol(X)) {
    X[,i] <- X[,i] - mean(X[,i])
    X[,i] <- X[,i]/sqrt(sum(X[,i]^{2}))
  }
  
  #  produce tilted correlation matrix 
  tilt <- diag(1, ncol(X))
  
  for (i in 1:ncol(X)){
    for (j in 1:ncol(X)){
      if (i != j){
        
        #  define new (x,y) variables 
        x <- X[,i]
        y <- X[,j]
        Xk <- X[,-c(i,j)]
        Z <- matrix(rep(1, nrow(X)), ncol = 1)
        
        #  find highly associated variables 
        for (k in 1:ncol(Xk)) {
          if ( abs(t(x) %*% Xk[,k]) > r 
               || abs(t(y) %*% Xk[,k]) > r) {
            Z <- cbind(Z, Xk[,k])
          }
        }
        
        if (ncol(Z)>1) {
          
          #  if highly associated varaibles exist 
          #  define projection matrix 
          if (det(solve(t(Z) %*% Z)) == 0) stop("choose a higher threshold (r), 
                                                the projection matrix is singular!")
          Z <- Z[,-1]
          P <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
          M <- diag(1,nrow(X)) - P
          
          #  transfrom x 
          x_str <- M %*% x
          a <- t(x) %*% M %*% x 
          
          #  transfrom y 
          y_str <- M %*% y
          b <- t(y) %*% M %*% y 
          
          #  replace off diagonals 
          tilt[i,j] = (sqrt(a*b))^{-1} * (t(x_str) %*% y_str)
          
        } else {
          
          #  if no highly associated vvaraibles exist 
          #  use straight correlation 
          tilt[i,j] = t(x) %*% y
          
        }
      }
    }
  }
  #  return tilted correlation
  return(tilt)
}