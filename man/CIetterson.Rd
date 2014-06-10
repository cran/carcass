\name{CIetterson}
\alias{CIetterson}
\title{
Confidence interval for the functions ettersonEq14, ettersonEq14v1 and ettersonEq14v2
}
\description{
Calculates Monte Carlo confidence intervals for the estimated carcass detection probability
when persistence probability and searcher efficiency are uncertain. 
}
\usage{
CIetterson(s, s.lwr, s.upr, f, f.lwr, f.upr, J, s.time.variance = "carcass age", 
f.time.variance = "number of searches", nsim = 1000, ci = 0.95)
}

\arguments{
  \item{s}{point estiate for persistence probability (see help file for functions etterson14, 
  ettersonEq14v1 or ettersonEq14v2)
}
  \item{s.lwr}{lower limit of the 95\% confidence interval of persistence probability
}
  \item{s.upr}{upper limit of the 95\% confidence interval of persistence probability
}
  \item{f}{point estimate for the searcher efficiency (see help file for functions etterson14, 
  ettersonEq14v1 or ettersonEq14v2
}
  \item{f.lwr}{lower limit of the 95\% confidence interval of searcher efficiency
}
  \item{f.upr}{upper limit of the 95\% confidence interval of searcher efficiency
}
  \item{J}{
vector of search intervals
}
  \item{s.time.variance}{
character, one of "date" or "carcass age"
}
  \item{f.time.variance}{
character, one of "date" or  "number of searches"
}
  \item{nsim}{
number of Monte Carlo simulations
}
  \item{ci}{
size of the confidence interval, default is 0.95
}
}
\details{
The time variance in s and f is either both with date or both with carcass age and 
number of searches, respectively.
In case of constant s and f, the function uses ettersonEq14 independent of the 
arguments s.time.variance or f.time.variance, when only one value is given for both 
parameters.
}
\value{
a list
 \item{p.lower}{lower limit of the confidence interval}
 \item{p.upper}{upper limit of the confidence interval}

}

\author{
F. Korner
}


\examples{
J <- c(2,3,2,4,3,5,3,2,3,4)
s <- plogis(seq(0.2, 2, length=sum(J)))
f <- plogis(seq(1.5, 0.9, length=length(J)))

s.lwr<- plogis(seq(0.2, 2, length=sum(J))-0.5)
f.lwr <- plogis(seq(1.5, 0.9, length=length(J))-0.3)

s.upr <- plogis(seq(0.2, 2, length=sum(J))+0.5)
f.upr <- plogis(seq(1.5, 0.9, length=length(J))+0.3)

CIetterson(s=s, s.lwr=s.lwr, s.upr=s.upr, f=f, f.lwr=f.lwr, f.upr=f.upr, J=J, nsim=100)
 # nsim is too low, please, increase!
}

\keyword{misc}
