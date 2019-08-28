# generate design matrix as VAR(1) process
# p = 500, n = 500
A = diag(0.4,1000)
S = toeplitz(GeometricSequence(1000,1,0.8))
X = as.matrix(VAR.sim(A, 500, include = 'none', varcov = S))

# generate y - the first five variables are relevant
e <- arima.sim(list(order=c(1,0,0), ar=.6), n=500)
y <- X[,200]-X[,400]+X[,600]-X[,800]+X[,1000] + e

# standardize simulated data 
y <- standardize(y)
X <- apply(X, 2, standardize)

# extract highly associated irrelevant varaibles 
to_keep <- which(rank(cor(X,y)) <= 55)
to_keep <- c(seq(200,1000,200),to_keep[!to_keep %in% seq(200,1000,200)][1:45])
print(to_keep)

# generate plot for straight correlations
# set simulation variables
sim_size <- seq(100,500,10) ; u <- length(sim_size)
p <- length(to_keep)
W <- matrix(rep(1,u*p), nrow = u, ncol = p)
for (i in 1:p){
  for (j in 1:length(sim_size)){
    # apply straight correlation for each sample size 
    y1 <- standardize(y[1:sim_size[j]])/sqrt(sim_size[j])
    x1 <- standardize(X[1:sim_size[j],to_keep[i]])/sqrt(sim_size[j])
    W[j,i] <- abs(cor(x1,y1))
  }
}

# plot straight correlations against sample size 
# colour relevant variables in red 
plot(W[,1]~sim_size, 
     type = 'l', 
     col = 'red', 
     ylim = c(0,max(W)), 
     xlab = 'n', 
     ylab = 'straight correlation', 
     main = 'p = 1000, n = 500')
for(i in 2:p){
  if (i %in% 2:5) {
    lines(W[,i]~sim_size, type = 'l', col = 'red')
  } else {
    lines(W[,i]~sim_size, type = 'l', col = 'grey')
  }
}

# generate threshold sequence
sim_size <- seq(100,500,10) # sequence of sub-samples
thr.seq <- c()
for (i in 1:length(sim_size)){
  thr.seq <- c(thr.seq, threshold(X[1:sim_size[i],]))
  cat('calculation', 100*round(i/length(sim_size),2), ' percent complete','\n')
}

# generate plot for tilted correlations
# set simulation variables
W <- matrix(rep(1,u*p), nrow = u, ncol = p)
for (i in 1:p){
  for (j in 1:length(sim_size)){
    # apply straight correlation for each sample size 
    y1 <- standardize(y[1:sim_size[j]])/sqrt(sim_size[j])
    X1 <- standardize(X[1:sim_size[j],])/sqrt(sim_size[j])
    W[j,i] <- abs(tilt(X1[,to_keep[i]],y1,X1[,-to_keep[i]], thr.seq[j] + stab(j)))
  }
  cat('calculations for variable', i, 'complete', '\n')
}

# plot tilted correlations against sample size 
# colour relevant variables in red 
plot(W[,1]~sim_size, 
     type = 'l', 
     col = 'red', 
     ylim = c(0,max(W)), 
     xlab = 'n', 
     ylab = 'tilted correlation', 
     main = 'p = 1000, n = 500')
for(i in 2:p){
  if (i %in% 2:5) {
    lines(W[,i]~sim_size, type = 'l', col = 'red')
  } else {
    lines(W[,i]~sim_size, type = 'l', col = 'grey')
  }
}