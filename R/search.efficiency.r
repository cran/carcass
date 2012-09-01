###########################################################################
# function to estimate searcher efficiency ith 95% confidence interval 
# based on experimental data
###########################################################################

# load libraries
library(arm)

search.efficiency <- function(dat, nsim=1000){
# dat data.frame containing the following variables: 
 #$ person     : names of the persons who searched
 #$ visibility : visibility class
 #$ detected   : number of detected items
 #$ notdetected: number of not detected items
 # nsim: number of simulations to be drawn from the posterior distributions to describe the 95% credible intervals
 
dat$visibility <- factor(dat$visibility)   
dat$person <- factor(dat$person)

npersons <- nlevels(dat$person)
nvisclass <- nlevels(dat$visibility)

if(npersons>2 & nvisclass>1){
  mod <- glmer(cbind(detected, notdetected) ~ visibility + (1|person), data=dat, family=binomial)
  newdat <- expand.grid(visibility=levels(dat$visibility), person=levels(dat$person))
  b <- fixef(mod)
  bperson <- matrix(b, nrow=npersons, ncol=3, byrow=TRUE) + cbind(ranef(mod)$person, rep(0, npersons), rep(0, npersons))
  newdat$f <- NA
  Xmat <- model.matrix(~visibility, data=newdat[newdat$person==levels(newdat$person)[1],]) 
  for(i in levels(dat$person)) newdat$f[newdat$person==i] <- plogis(Xmat%*%as.numeric(bperson[i,]))
   
  bsim <- sim(mod, n.sim=nsim)
  
  predmat <- matrix(nrow=nrow(newdat), ncol=nsim)
  for(i in levels(dat$person)){ 
    for(j in 1:nsim) {
    bpersi <- b + c(bsim@ranef$person[j,i,1],0, 0)
    predmat[newdat$person==i, j] <- plogis(Xmat%*%bpersi)
    }
    }
  newdat$lwr <- apply(predmat, 1, quantile, prob=0.025)
  newdat$upr <- apply(predmat, 1, quantile, prob=0.975)
  
  overall <- data.frame(visibility=levels(dat$visibility))
  Xmat <- model.matrix(~visibility, data=overall)
  overall$f <- plogis(Xmat%*%b)
  predmat <- matrix(nrow=nrow(overall), ncol=nsim)
  for(i in 1:nsim) predmat[,i] <- plogis(Xmat%*%bsim@fixef[i,])
  overall$lwr <- apply(predmat, 1, quantile, prob=0.025)
  overall$upr <- apply(predmat, 1, quantile, prob=0.975)
  }



if(npersons==2 & nvisclass>1){
  mod <- glmer(cbind(detected, notdetected) ~ visibility + person, data=dat, family=binomial)
  newdat <- expand.grid(visibility=levels(dat$visibility), person=levels(dat$person))
  b <- fixef(mod)
  Xmat <- model.matrix(~visibility + person, data=newdat) 
  newdat$f <- plogis(Xmat%*%b)
   
  bsim <- sim(mod, n.sim=nsim)
  
  predmat <- matrix(nrow=nrow(newdat), ncol=nsim)
  for(i in 1:nsim) predmat[, i] <- plogis(Xmat%*%bsim@fixef[i,])

  newdat$lwr <- apply(predmat, 1, quantile, prob=0.025)
  newdat$upr <- apply(predmat, 1, quantile, prob=0.975)
  
  overall <- data.frame(visibility=levels(dat$visibility))
  predmat1 <- matrix(nrow=nrow(overall), ncol=nsim)
  for(i in 1:nsim) predmat1[,i] <- tapply(predmat[,i], newdat$visibility, mean)
  overall$f <- apply(predmat1, 1, mean)
  overall$lwr <- apply(predmat1, 1, quantile, prob=0.025)
  overall$upr <- apply(predmat1, 1, quantile, prob=0.975)
  }
     

if(npersons==1 & nvisclass>1){
  mod <- glmer(cbind(detected, notdetected) ~ visibility, data=dat, family=binomial)
  newdat <- expand.grid(visibility=levels(dat$visibility), person=levels(dat$person))
  b <- fixef(mod)
  Xmat <- model.matrix(~visibility, data=newdat) 
  newdat$f <- plogis(Xmat%*%b)
   
  bsim <- sim(mod, n.sim=nsim)
  
  predmat <- matrix(nrow=nrow(newdat), ncol=nsim)
  for(i in 1:nsim) predmat[, i] <- plogis(Xmat%*%bsim@fixef[i,])

  newdat$lwr <- apply(predmat, 1, quantile, prob=0.025)
  newdat$upr <- apply(predmat, 1, quantile, prob=0.975)
  
  overall <- newdat
  }




if(npersons>2 & nvisclass==1){
  mod <- glmer(cbind(detected, notdetected) ~ 1 + (1|person), data=dat, family=binomial)
  newdat <- expand.grid(visibility=levels(dat$visibility), person=levels(dat$person))
  b <- fixef(mod)
  bperson <- matrix(b, nrow=npersons, ncol=3, byrow=TRUE) + cbind(ranef(mod)$person, rep(0, npersons), rep(0, npersons))
  newdat$f <- NA
  Xmat <- model.matrix(~1, data=newdat[newdat$person==levels(newdat$person)[1],]) 
  for(i in levels(dat$person)) newdat$f[newdat$person==i] <- plogis(Xmat%*%as.numeric(bperson[i,]))
   
  bsim <- sim(mod, n.sim=nsim)
  
  predmat <- matrix(nrow=nrow(newdat), ncol=nsim)
  for(i in levels(dat$person)){ 
    for(j in 1:nsim) {
    bpersi <- b + c(bsim@ranef$person[j,i,1],0, 0)
    predmat[newdat$person==i, j] <- plogis(Xmat%*%bpersi)
    }
    }
  newdat$lwr <- apply(predmat, 1, quantile, prob=0.025)
  newdat$upr <- apply(predmat, 1, quantile, prob=0.975)
  
  overall <- data.frame(visibility=levels(dat$visibility))
  Xmat <- model.matrix(~1, data=overall)
  overall$f <- plogis(Xmat%*%b)
  predmat <- matrix(nrow=nrow(overall), ncol=nsim)
  for(i in 1:nsim) predmat[,i] <- plogis(Xmat%*%bsim@fixef[i,])
  overall$lwr <- apply(predmat, 1, quantile, prob=0.025)
  overall$upr <- apply(predmat, 1, quantile, prob=0.975)
  }



if(npersons==2 & nvisclass==1){
  mod <- glmer(cbind(detected, notdetected) ~ 1 + person, data=dat, family=binomial)
  newdat <- expand.grid(visibility=levels(dat$visibility), person=levels(dat$person))
  b <- fixef(mod)
  Xmat <- model.matrix(~1 + person, data=newdat) 
  newdat$f <- plogis(Xmat%*%b)
   
  bsim <- sim(mod, n.sim=nsim)
  
  predmat <- matrix(nrow=nrow(newdat), ncol=nsim)
  for(i in 1:nsim) predmat[, i] <- plogis(Xmat%*%bsim@fixef[i,])

  newdat$lwr <- apply(predmat, 1, quantile, prob=0.025)
  newdat$upr <- apply(predmat, 1, quantile, prob=0.975)
  
  overall <- data.frame(visibility=levels(dat$visibility))
  predmat1 <- matrix(nrow=nrow(overall), ncol=nsim)
  for(i in 1:nsim) predmat1[,i] <- tapply(predmat[,i], newdat$visibility, mean)
  overall$f <- apply(predmat1, 1, mean)
  overall$lwr <- apply(predmat1, 1, quantile, prob=0.025)
  overall$upr <- apply(predmat1, 1, quantile, prob=0.975)
  }
     

if(npersons==1 & nvisclass==1){
  mod <- glmer(cbind(detected, notdetected) ~ 1, data=dat, family=binomial)
  newdat <- expand.grid(visibility=levels(dat$visibility), person=levels(dat$person))
  b <- fixef(mod)
  Xmat <- model.matrix(~1, data=newdat) 
  newdat$f <- plogis(Xmat%*%b)
   
  bsim <- sim(mod, n.sim=nsim)
  
  predmat <- matrix(nrow=nrow(newdat), ncol=nsim)
  for(i in 1:nsim) predmat[, i] <- plogis(Xmat%*%bsim@fixef[i,])

  newdat$lwr <- apply(predmat, 1, quantile, prob=0.025)
  newdat$upr <- apply(predmat, 1, quantile, prob=0.975)
  
  overall <- newdat
  }

  return(list(f.perperson=newdat, f.average=overall))
  }
  



