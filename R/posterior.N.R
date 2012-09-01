posterior.N <-
function(p, nf=0, maxN=50, ci.int=0.95, plot=TRUE, dist=FALSE){
# p = probability that a killed animal is detected by a seacher
# nf = number of carcasses found
# maxN = maximal possible number of fatalities
# ci.int = size of the credible interval that should be calculated
# plot: posterior distribution is plotted if TRUE 
# dist: posterior distribution is given if  TRUE
#---------------------------------------------------------------
N <- nf:maxN
if(nf==0) pN <- p*(1-p)^(N-nf)
if(nf>0) {
  denom <- sum(choose(N, nf) * (1-p)^(N-nf))
  pN <- choose(N, nf)*(1-p)^(N-nf)/denom
  pN <- c(rep(0, nf), pN)
  N <- c(rep(0, nf), N)  
}

if(plot) plot(N, pN, type="h", lwd=5, lend="butt", xlab="Number of fatalities", ylab="Posterior density")
index <- cumsum(pN)<ci.int
indexLower <- cumsum(pN)<(1-ci.int)/2
indexUpper <- cumsum(pN)<1-(1-ci.int)/2
if(nf==0) interval <- c(nf, min(N[!index]))   
if(nf>0)  interval <- c(min(N[!indexLower]), min(N[!indexUpper])) 
if(interval[2]==Inf) cat("Upper limit of CI larger than maxN! -> increase maxN\n")
expected <- min(N[!cumsum(pN)<0.5])
results <- list(interval=interval, expected=expected)
if(dist==TRUE) results <- list(interval=interval, expected=expected, pN=pN)
results
}


