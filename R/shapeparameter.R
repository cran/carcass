shapeparameter <-
function(m, lwr=NA, upr=NA, se=NA){
# m = estimate
# lwr, upr = lower and upper limit of the 95 % credible or confidence 
# se = standard error of the proportion
# interval
#--------------------------------------------------------------------
if(is.na(se)){
  ci <- upr - lwr
  sigma2 <-(ci/4)^2
  }
if(!is.na(se)){
  sigma2 <- se^2
  }
a <- m*(m*(1-m)/sigma2-1)
b <- (1-m)*(m*(1-m)/sigma2-1)
return(list(a=a,b=b))
}

