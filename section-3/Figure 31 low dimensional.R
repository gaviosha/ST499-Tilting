# generate design matrix as VAR(1) process
# p = 50, n = 500
A = diag(0.4,50)
S = toeplitz(GeometricSequence(50,1,0.8))
X = as.matrix(VAR.sim(A, 500, include = 'none', varcov = S))

# generate y - the first five variables are relevant
e <- arima.sim(list(order=c(1,0,0), ar=.6), n=500) 
y <- X[,10]-X[,20]+X[,30]-X[,40]+X[,50] + e

# standardize simulated data 
y <- standardize(y)
X <- apply(X, 2, standardize)

# generate threshold sequence
sim_size <- seq(100,500,10) # sequence of sub-samples
thr.seq <- c()
for (i in 1:length(sim_size)){
  thr.seq <- c(thr.seq, threshold(X[1:sim_size[i],]))
  cat(round(i/length(sim_size),1), ' complete','\n')
}

# generate plot for tilted correlations 
# set simulation variables
u <- length(sim_size)
p <- 50 
W <- matrix(rep(1,u*p), nrow = u, ncol = p)

for (i in 1:p){
  for (j in 1:length(sim_size)){
    # apply tilted correlation for each sample size 
    y1 <- standardize(y[1:sim_size[j]])/sqrt(sim_size[j])
    X1 <- apply(X[1:sim_size[j], ],2,standardize)/sqrt(sim_size[j])
    W[j,i] <- abs(tilt(X1[,i],y1,X1[,-i], thr.seq[j]))
  }
}

# plot tilted correlations against sample size 
# colour relevant variables in red 
plot(W[,1]~sim_size, 
     type = 'l', 
     col = 'grey', 
     ylim = c(0,max(W)), 
     xlab = 'n', 
     ylab = 'tilted correlation', 
     main = 'p = 50, n = 500')
for(i in 2:p){
  if (i%%10 == 0) {
    lines(W[,i]~sim_size, type = 'l', col = 'red')
  } else {
    lines(W[,i]~sim_size, type = 'l', col = 'grey')
  }
}

# generate plot for straight correlations 
# set simulation variables
W <- matrix(rep(1,u*p), nrow = u, ncol = p)

for (i in 1:p){
  for (j in 1:length(sim_size)){
    # apply straight correlation for each sample size 
    y1 <- standardize(y[1:sim_size[j]])/sqrt(sim_size[j])
    x1 <- standardize(X[1:sim_size[j],i])/sqrt(sim_size[j])
    W[j,i] <- abs(cor(x1,y1))
  }
}

# plot straight correlations against sample size 
# colour relevant variables in red 
plot(W[,1]~sim_size, 
     type = 'l', 
     col = 'grey', 
     ylim = c(0,max(W)), 
     xlab = 'n', 
     ylab = 'straight correlation', 
     main = 'p = 50, n = 500')
for(i in 2:p){
  if (i%%10 == 0) {
    lines(W[,i]~sim_size, type = 'l', col = 'red')
  } else {
    lines(W[,i]~sim_size, type = 'l', col = 'grey')
  }
}