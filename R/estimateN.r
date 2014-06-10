# Function to estimate the number of killed animals from carcass counts 
# when the detection probability is estimated with uncertainty (95% CI)
# Author: Fraenzi Korner-Nievergelt
# This function wraps the code given in Appendix II of the following paper where you find detailed comments in each step:

# Korner-Nievergelt F, Korner-Nievergelt P, Behr O, Niermann I, Brinkmann R & Hellriegel B (2011)
# A new method to determine bird and bat fatality at wind energy turbines from carcass searches. Wildl. Biol. 17: 1-14
#-------------------------------------------------------------------------------



# R-Code written and tested for R. 2.12.0
#-------------------------------------------------------------------------------

#Estimating the number of fatalities with a credible interval
#Step 1:
# Describe the uncertainty in the estimates for searcher efficiency and carcass persistence probability by a beta-distribution, i.e. transform the 95 % CI into the shape parameters of the beta-distributions.

estimateN <- function(count, p=NA, p.lower=NA, p.upper=NA, f=NA, f.lower=NA, f.upper=NA, 
                        s=NA, s.lower=NA, s.upper=NA, a=1, a.lower=1, a.upper=1,
                        pform="korner", d=1, n=NA, J=NA, maxn=1000, nsim=1000, plot=TRUE, postdist=FALSE, k=1){

# count   number of carcasses found
# p       estimate for detection probability
# p.lower lower limit of 95% CI of detection probability
# p.upper upper limit of 95% CI of detection probability

# eihter p or f and s have to be provided!!

# f       estimate for searcher efficiency
# f.lower lower limit of 95% CI of searcher efficiency
# f.upper upper limit of 95% CI of searcher efficiency
# s       estimate for persistence probability
# s.lower lower limit of 95% CI of persistence probability
# s.upper upper limit of 95% CI of persistence probability
# a       estimate for the proportion of killed animals that falled into a searched area
# a.lower lower limit of 95% CI of proportion of killed animals in searched area
# a.upper upper limit of 95% CI of proportion of killed animals in searched area
# d       search interval, the number days between two searches
# n       number of searches
# J       vector with the lengths of the search intervals (if search intervals are irregular)
# pform   formula used to estimate p, one of "korner", "huso", "erickson", "etterson"

# maxn   maximal possible number of animals killed for which the posterior probability is estimated (should not be too high but also not be too low!)
# nsim   number of Monte Carlo simulations


if(is.finite(f)[1]){
f.a <- shapeparameter(f, f.lower, f.upper)$a
f.b <- shapeparameter(f, f.lower, f.upper)$b
s.a <- shapeparameter(s, s.lower, s.upper)$a
s.b <- shapeparameter(s, s.lower, s.upper)$b
if(a.lower<a.upper){
a.a <- shapeparameter(a, a.lower, a.upper)$a
a.b <- shapeparameter(a, a.lower, a.upper)$b
}

Npostdist <- numeric(maxn+1) 

for(i in 1:nsim){

fr <- rbeta(length(f), f.a, f.b)
sr <- rbeta(length(s), s.a, s.b)
ar <- ifelse(a.lower<a.upper, rbeta(1, a.a, a.b), a)

if(pform=="korner") pr <- pkorner(s=sr, f=fr, d=d, n=n, k=k, search.efficiency.constant=ifelse(k==1, TRUE, FALSE))
if(pform=="huso") pr <- phuso(s=sr, f=fr, d=d)
if(pform=="erickson") pr <- perickson(t.bar=-1/log(sr), f=fr, d=d)  # 
if(pform=="etterson") pr <- ettersonEq14v2(s=sr, f=fr, J=J)

postNtemp <- posteriorN(nf=count, p=pr*ar, maxN=maxn,  plot=FALSE, dist=TRUE)

Npostdist <- Npostdist + postNtemp$pN    #Sum the posterior densities over all nsim simulations.

} # close loop i
if(pform=="korner") pm <- pkorner(s=s, f=f, d=d, n=n, k=k, search.efficiency.constant=ifelse(k==1, TRUE, FALSE))
if(pform=="huso") pm <- phuso(s=s, f=f, d=d)
if(pform=="erickson") pm <- perickson(t.bar=-1/log(s), f=f, d=d)  # 
if(pform=="etterson") pm <- ettersonEq14v2(s=s, f=f, J=J)
HT.estimate <- count/(pm*a) # Horwitz-Thompson estimate
} # close sim based on f and s


if(is.finite(p)){
p.a <- shapeparameter(p, p.lower, p.upper)$a
p.b <- shapeparameter(p, p.lower, p.upper)$b
HT.estimate <- count/(p*a) # Horwitz-Thompson estimate

Npostdist <- numeric(maxn+1) 


for(i in 1:nsim){
pr <- rbeta(1, p.a, p.b)

postNtemp <- posteriorN(nf=count, p=pr, maxN=maxn,  plot=FALSE, dist=TRUE)
Npostdist <- Npostdist + postNtemp$pN    #Sum the posterior densities over all nsim simulations.

} # close loop i
} # close sim based on p


Npostdist.sc <- Npostdist/nsim
indexLower <- cumsum(Npostdist.sc) < 0.025
indexMedian <- cumsum(Npostdist.sc) < 0.5 
indexUpper <- cumsum(Npostdist.sc) < 0.975
lower <- min(c(0:maxn)[!indexLower])  
estimate.median <- min(c(0:maxn)[!indexMedian])  
#estimate.mean <- sum(c(0:maxn)*Npostdist.sc)
upper <- min(c(0:maxn)[!indexUpper])

if(plot) plot(0:maxn, Npostdist.sc , type="h", lwd=5, lend="butt", xlab="Number of fatalities", ylab="Posterior density", xlim=c(0,50))
if(!postdist) {
result <- c(estimate.median, lower, upper, HT.estimate)
names(result) <- c("estimate", "lower", "upper", "HT.estimate")
}
if(postdist) {
result <- list(estimate=estimate.median, lower=lower, upper=upper, HT.estimate=HT.estimate, 
               postdist=Npostdist.sc)
}

return(result)
}






