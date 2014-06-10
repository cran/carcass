# Log of changes:
# - FK, 23.06.2013: Bug when using search.efficiency.constant=FALSE removed.
# - TR, 26.06.2013: Function can also handle non-constant persistence probability.
# - FK, 26.08.2013: fixed the problem that function returned NA when for non -constant persistence probabilities less number of estimates than n*d were available. The function now assumes that the proportion of persistant carcasses is zero after the last estimate. 
# - FK, 26.10.2013: inserted for non-constant persistence prob the line for the first search (at day =1)

pkorner <-
function(s, s.lower=NA, s.upper=NA, f, f.lower=NA, f.upper=NA, d, n, k=0.25, search.efficiency.constant=TRUE, CI=FALSE, nsim=1000){
  # s = probability that a carcass remains 24 hours (scalar if constant, vector otherwise)
  # f = probability that a carcass is detected by a searcher during a search given it persisted to the search
  # d = (average) number of days between two searches
  # n = number of searches (n * d = length of study period)
  # k = factor by which the search efficiency is multiblied after each search
  #--------------------------------------------------------------
  persistence.probability.constant <- length(s)==1
  if(search.efficiency.constant & persistence.probability.constant){
    x <- (1-f)*s^d
    A <- s*(1-s^d)/(1-s)
    summep <- numeric(n)
    for(i in 0:(n-1)) summep[i+1] <- (n-i)*x^i
    p <- A*f*sum(summep)/(d*n)
  }
  if(!search.efficiency.constant & persistence.probability.constant){
    A<-s*(1-s^d)/(1-s)
    summepfound<-A*f
    for(x in 2:n){
      summep<-1
      for(j in 1:(x-1)) summep<-summep+k^(x-j) * s^((x-j)*d) * productsef(f, k, n=x, j)
      summepfound<-summepfound+A*f*summep
    }
    p<-summepfound/(d*n)     
  }
  if(!persistence.probability.constant){
    f <- rep(f, n)
    if(!search.efficiency.constant) for(j in 2:length(f)) f[j] <- f[j-1]*k
    c <- array(0, dim=c(n*d, n*d))
    found <- array(0, dim=c(n*d, n*d))
    diag(c) <- s[1]
    ns <- length(s)      # fk inserted
    if(ns<10) message(paste("Proportion of persisting carcasses is only given for", ns, "days. It is assumed that no carcass persists longer than", ns, "days!"))
    if(d==1) found[1,1] <- f[1]*c[1,1]
    for(day in 2:(n*d)) {
      for(cohort in 1:(day-1)) c[day,cohort] <- (c[day-1,cohort] - found[day-1,cohort])*ifelse((day-cohort+1)<= ns, s[day-cohort+1], 0) # fk inserted ifelse function
      for(cohort in 1:day) if(day%%d==0) found[day,cohort] <- c[day,cohort]*f[ceiling((day-cohort+1)/d)]
    }
    p <- sum(found)/(n*d)
  }
  
### Get confidence intervals
  if(CI){
    if(is.na(s.lower[1]) | is.na(s.upper[1]) | is.na(f.lower) | is.na(f.upper)) stop("'Please provide CI for s and f.")
    f.a <- shapeparameter(f, f.lower, f.upper)$a
    f.b <- shapeparameter(f, f.lower, f.upper)$b
    s.a <- shapeparameter(s, s.lower, s.upper)$a
    s.b <- shapeparameter(s, s.lower, s.upper)$b
    resCI <- numeric(nsim)
    for(i in 1:nsim) {
      resCI[i] <- pkorner(s=rbeta(length(s.a), s.a, s.b), f=rbeta(1,f.a, f.b), d=d, n=n, search.efficiency.constant=search.efficiency.constant, CI=FALSE, nsim=nsim)
    }
    names(p) <- "p"
    p <- c(p, quantile(x=resCI, probs=c(0.025, 0.975)))
  }
  
### Return results  
  p
}

