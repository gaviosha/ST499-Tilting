# simulates average tilted correlation
# using method 3.1
# requires package ppcor 

require(ppcor)

intersection_sim_mean <- function(sims, n, p, Rmin, Rmax) {
  
  for (i in 1:sims) {
    
    cat("Beginning intersection simulation:", i)
    cat("\n")
    
    #  simulate matrices 
    W <- rnorm(n)
    for (j in 1:(p-1)) {
      W <- cbind(W, rnorm(n))
    }
    
    #  threshold and holder 
    R <- seq(from = Rmin, to = Rmax, by = 0.01)
    xx <- c()
    
    for (k in 1:length(R)) {
      xx <- c(xx, mean(abs(intersection_tilt(W, R[k]) - diag(1,p))))
      cat(round(k/length(k), digits = 2), "% complete")
      cat("\n")
    }
    
    #  plot results 
    plot(R,xx, type = "l", 
         xlab = "Threshold", 
         ylab = "maximum correlation")
    lines(R, rep(mean(abs(cor(W) - diag(1,p))), length(R)), col = "red")
    lines(R, rep(mean(abs(pcor(W)$estimate - diag(1,p))), length(R)), col = "green")
    
    cat("\n")
    cat("\n")
  }
}