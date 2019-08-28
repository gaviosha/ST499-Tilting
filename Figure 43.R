###################################
# Heavy tails analysis : R2 = 0.9 #
###################################

set.seed(123) # set seed 

# generate design matrix as VAR(1) process
# p = 1000, n = 200
A = diag(0.4,1000)
S = toeplitz(GeometricSequence(1000,1,0.8))
X = as.matrix(VAR.sim(A, 200, include = 'none', 
                      innov = rmvt(n = 200, sigma = S, df = 3)))

# work out sigma for R2 = 0.9
b = rep(0,1000)
b[c(200,400,600,800,1000)] = c(1,-1,1,-1,1)
s = (t(b)%*%S%*%b/(1-0.4^2))*((1/0.9)-1)*(1-0.6^2)

# generate y - the first five variables are relevant
v =rnorm(n = 200, sd = sqrt(s))
e <- arima.sim(list(order=c(1,0,0), ar=.6), n=200, innov = v)
y <- X[,200]-X[,400]+X[,600]-X[,800]+X[,1000] + e

# generate threshold once 
start <- Sys.time()
thr.for.all <- threshold(X)
end <- Sys.time()
end - start 
print(thr.for.all)

print('Now starting: REGULAR TILTING')
active.tilt <- tilt.it.normal(X,y,20, thr = thr.for.all)
print(active.tilt)

print('Now starting: EFFICIENT TILTING')
active.efficient <- tilt.it.efficient(X,y,20, thr = thr.for.all)
print(active.efficient)

print('Now starting: STABILITY TILTING')
active.stable <- tilt.it.small.n.stable(X,y,20, thr = thr.for.all)
print(active.stable)

## generate FPR and TPR for all methods 
# tilting 
TPR.tilt <- make.TPR(active.tilt)
FPR.tilt <- make.FPR(active.tilt)
# efficient tilting 
TPR.efficient <- make.TPR(active.efficient)
FPR.efficient <- make.FPR(active.efficient)
# stable tilting 
TPR.stable <- make.TPR(active.stable)
FPR.stable <- make.FPR(active.stable)

## plot results 
# tilting in red 
plot(FPR.tilt,TPR.tilt, type = 'b', pch = 1,
     col = 'red', 
     xlab = 'FPR',
     xlim = c(0,max(c(FPR.tilt,FPR.efficient,FPR.stable))),
     ylab = 'TPR',
     ylim = c(0,max(c(TPR.tilt,TPR.efficient,TPR.stable)))
)
# efficient tilting in blue 
lines(FPR.efficient, TPR.efficient, type = 'b', pch = 2, col = 'blue')
# stability tilting in green 
lines(FPR.stable, TPR.stable, type = 'b', pch = 3, col = 'green')

# add legend 
# legend(0.01, 0.9, legend=c("Tilting", "Efficient Tilting","Stable Tilting"),
col=c("red", "blue","green"), lty=1,pch = 1:3)


###################################
# Heavy tails analysis : R2 = 0.6 #
###################################

set.seed(123) # set seed 

# generate design matrix as VAR(1) process
# p = 1000, n = 200
A = diag(0.4,1000)
S = toeplitz(GeometricSequence(1000,1,0.8))
X = as.matrix(VAR.sim(A, 200, include = 'none', 
                      innov = rmvt(n = 200, sigma = S, df = 3)))

# work out sigma for R2 = 0.6
b = rep(0,1000)
b[c(200,400,600,800,1000)] = c(1,-1,1,-1,1)
s = (t(b)%*%S%*%b/(1-0.4^2))*((1/0.6)-1)*(1-0.6^2)

# generate y - the first five variables are relevant
v =rnorm(n = 200, sd = sqrt(s))
e <- arima.sim(list(order=c(1,0,0), ar=.6), n=200, innov = v)
y <- X[,200]-X[,400]+X[,600]-X[,800]+X[,1000] + e

# seed is same, we can re-use threshold! 
print(thr.for.all)

print('Now starting: REGULAR TILTING')
active.tilt <- tilt.it.normal(X,y,20, thr = thr.for.all)
print(active.tilt)

print('Now starting: EFFICIENT TILTING')
active.efficient <- tilt.it.efficient(X,y,20, thr = thr.for.all)
print(active.efficient)

print('Now starting: STABILITY TILTING')
active.stable <- tilt.it.small.n.stable(X,y,20, thr = thr.for.all)
print(active.stable)

## generate FPR and TPR for all methods 
# tilting 
TPR.tilt <- make.TPR(active.tilt)
FPR.tilt <- make.FPR(active.tilt)
# efficient tilting 
TPR.efficient <- make.TPR(active.efficient)
FPR.efficient <- make.FPR(active.efficient)
# stable tilting 
TPR.stable <- make.TPR(active.stable)
FPR.stable <- make.FPR(active.stable)

## plot results 
# tilting in red 
plot(FPR.tilt,TPR.tilt, type = 'b', pch = 1,
     col = 'red', 
     xlab = 'FPR',
     xlim = c(0,max(c(FPR.tilt,FPR.efficient,FPR.stable))),
     ylab = 'TPR',
     ylim = c(0,max(c(TPR.tilt,TPR.efficient,TPR.stable)))
)
# efficient tilting in blue 
lines(FPR.efficient, TPR.efficient, type = 'b', pch = 2, col = 'blue')
# stability tilting in green 
lines(FPR.stable, TPR.stable, type = 'b', pch = 3, col = 'green')


###################################
# Heavy tails analysis : R2 = 0.3 #
###################################

set.seed(123) # set seed 

# generate design matrix as VAR(1) process
# p = 1000, n = 200
A = diag(0.4,1000)
S = toeplitz(GeometricSequence(1000,1,0.8))
X = as.matrix(VAR.sim(A, 200, include = 'none', 
                      innov = rmvt(n = 200, sigma = S, df = 3)))

# work out sigma for R2 = 0.3
b = rep(0,1000)
b[c(200,400,600,800,1000)] = c(1,-1,1,-1,1)
s = (t(b)%*%S%*%b/(1-0.4^2))*((1/0.3)-1)*(1-0.6^2)

# generate y - the first five variables are relevant
v =rnorm(n = 200, sd = sqrt(s))
e <- arima.sim(list(order=c(1,0,0), ar=.6), n=200, innov = v)
y <- X[,200]-X[,400]+X[,600]-X[,800]+X[,1000] + e

# seed is same, we can re-use threshold! 
print(thr.for.all)

print('Now starting: REGULAR TILTING')
active.tilt <- tilt.it.normal(X,y,20, thr = thr.for.all)
print(active.tilt)

print('Now starting: EFFICIENT TILTING')
active.efficient <- tilt.it.efficient(X,y,20, thr = thr.for.all)
print(active.efficient)

print('Now starting: STABILITY TILTING')
active.stable <- tilt.it.small.n.stable(X,y,20, thr = thr.for.all)
print(active.stable)

## generate FPR and TPR for all methods 
# tilting 
TPR.tilt <- make.TPR(active.tilt)
FPR.tilt <- make.FPR(active.tilt)
# efficient tilting 
TPR.efficient <- make.TPR(active.efficient)
FPR.efficient <- make.FPR(active.efficient)
# stable tilting 
TPR.stable <- make.TPR(active.stable)
FPR.stable <- make.FPR(active.stable)

## plot results 
# tilting in red 
plot(FPR.tilt,TPR.tilt, type = 'b', pch = 1,
     col = 'red', 
     xlab = 'FPR',
     xlim = c(0,max(c(FPR.tilt,FPR.efficient,FPR.stable))),
     ylab = 'TPR',
     ylim = c(0,max(c(TPR.tilt,TPR.efficient,TPR.stable)))
)
# efficient tilting in blue 
lines(FPR.efficient, TPR.efficient, type = 'b', pch = 2, col = 'blue')
# stability tilting in green 
lines(FPR.stable, TPR.stable, type = 'b', pch = 3, col = 'green')