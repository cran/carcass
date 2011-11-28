shapeparameter <-
function(m, lwr, upr){
# m = estimate
# lwr, upr = lower and upper limit of the 95 % credible or confidence 
# interval
#--------------------------------------------------------------------
ci <- upr - lwr
sigma2 <-(ci/4)^2
a <- m*(m*(1-m)/sigma2-1)
b <- (1-m)*(m*(1-m)/sigma2-1)
list(a=a,b=b)
}

