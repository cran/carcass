\name{estimateN}
\alias{estimateN}
\title{Estimation of number of killed animals based on carcass searches and estimates 
for detection probability
}
\description{
Estimates the number of killed animals from carcass searches taking into account 
the uncertainty in the estimates for carcass detection probability (or alternatively 
for searcher efficiency and carcass persistence probability).
}
\usage{
estimateN(count, p = NA, p.lower = NA, p.upper = NA,
          f = NA, f.lower = NA, f.upper = NA, s = NA, 
          s.lower = NA, s.upper = NA, arrival="discrete", a=1, a.lower=1, a.upper=1,
          pform = "korner", d = 1, n = NA, J=NA, maxn = 1000, 
          nsim = 1000, plot = TRUE, postdist=FALSE, k=1,  x=c(1:10))
}
\arguments{
  \item{count}{
number of carcasses found
}
  \item{p}{
estimate for detection probability (alternatively give f, s, d and n)
}
  \item{p.lower}{
lower limit of 95\% CI of detection probability
}
  \item{p.upper}{
upper limit of 95\% CI of detection probability
}
  \item{f}{
estimate for searcher efficiency
}
  \item{f.lower}{
lower limit of 95\% CI of searcher efficiency
}
  \item{f.upper}{
upper limit of 95\% CI of searcher efficiency
}
  \item{s}{
estimate for persistence probability: this is a scalar in case of constant persistence probability.
In case of non-constant persistence probability, it is a vector. ATTENTION: 
If 'pform="korner"', the proportion
of remaining carcasses until each day after death (survivor function) have to be given. 
If 'pform="etterson"', the daily persistence probabilities for each day after death are needed.
}
  \item{s.lower}{
lower limit of 95\% CI of persistence probability
}
  \item{s.upper}{
upper limit of 95\% CI of persistence probability
}
  \item{arrival}{
either "uniform" or "discrete".
assumption of the distribution of carcass arrival times: 
if "uniform" it is assumed that carcasses arrive and are removed at a constant rate 
over time. 
if "discrete" (default) it is assumed that all carcasses arrive every day at the same time and 
and persistence probability given for the first day corresponds to the time from
death to the time when carcass searches normally take place on the first day (the search 
does not necessary need to take place on the first day).
}
  \item{a}{
estimate for the proportion of killed animals that falled into the searched area. 
default is 1 (complete area searched).
}
  \item{a.lower}{
lower limit of 95\% CI of the proportion of killed animals that falled into the 
searched area. default is 1 (complete area searched).
}
  \item{a.upper}{
upper limit of 95\% CI of the proportion of killed animals that falled into the 
searched area. default is 1 (complete area searched).
}

  \item{pform}{
formula used to estimate p, one of "korner", "huso", "erickson", "etterson"
}
  \item{d}{
 search interval, the number days (or any other time unit) between two searches
}
  \item{n}{
number of searches
}
  \item{J}{
vector of length of search intervals. This argument can only be used with pform="etterson". 
For the other methods regular search intervals are used and their length is given in the argument d.
}
  \item{maxn}{
maximal possible number of animals killed for which the posterior probability is estimated (should not be too high but also not be too low!)
}
  \item{nsim}{
number of Monte Carlo simulations
}
  \item{plot}{
logical, TRUE = posterior probability is plotted
}
  \item{postdist}{
logical, TRUE = the posterior probability distribution of the number of fatalities is given in the output. See details.
}
  \item{k}{
numeric, factor of decrease of searcher efficiency with number of searches, see help to pkorner.
}
  \item{x}{
vector of integers. The posterior probability of the number of killed animals being larger than x will be calculated.
}

}
\details{
Either p (with p.lower and p.upper) OR f (with f.lower and f.upper), 
s (with s.lower and s.upper), n and d are needed. In the latter case, p is calculated
from f, s, n and d using the formula specified in the argument 'pform'.

A documented code of this function if given in the Appendix of Korner-Nievergelt et al. (2011) J. Wildlife. Biol.

If postdist is TRUE, the output is a list, if postdist is FALSE, the output is a named numeric vector. 
}
\value{
If postdist=FALSE, the function produces a vector with the median of the posterior distribution ("estimate"),
the 2.5\% and the 97.5\% quantiles ("lower", and "upper") and the Horvitz-Thompson estimate ("HT.estimate", which equals
count/p). 
If postdist=TRUE, the output is a list with the following elements:
 \item{estimate}{median of the posterior distribution}
 \item{lower}{2.5\% quantile of the posterior distribution}
 \item{upper}{97.5\% quantile of the posterior distribution}
 \item{HT.estimate}{Horvitz-Thompson estimate, i.e. count/p}
 \item{P.true.larger.x}{Posterior probability that the true mortality is larger than x.}

}
\references{
Erickson WP, Jeffrey J, Kronner K, Bay K (2004) Stateline Wind Project Wildlife. 
Monitoring Final Report, July 2001 - December 2003. - Technical report peer-reviewed 
by and submitted to FPL Energy, the Oregon Energy Facility Siting Council, and the Stateline 
Technical Advisory Committee, Oregon, USA, 105 pp\cr
Etterson, MA (2013) Hidden Markov models for estimating animal mortality from antropogenic 
hazards. Ecological Applications 23: 1915-1925\cr
Huso M (2010) An estimator of wildlife fatality from observed carcasses. Environmetrics 
22: 318-329\cr
Korner-Nievergelt F, Korner-Nievergelt P, Behr O, Niermann I, Brinkmann R, Hellriegel B (2011) 
A new method to determine bird and bat fatality at wind energy turbines from carcass searches.
Wildlife Biology 17: 350-363
Huso M, Dathorp D in prep.
}
\author{
F. Korner-Nievergelt
}


\seealso{
\code{\link{posteriorN}}
}
\examples{

estimateN(count=3, f=0.72, f.lower=0.62, f.upper=0.81, s=0.84, s.lower=0.64, 
            s.upper=0.94, d=2, pform="korner", n=100, maxn=500, nsim=1000, 
            plot=TRUE)

estimateN(count=3, f=0.72, f.lower=0.62, f.upper=0.81, s=0.84, s.lower=0.64, 
            s.upper=0.94, d=2, pform="huso", maxn=500, nsim=1000, plot=TRUE)

res.p <- pkorner(f=0.72, f.lower=0.62, f.upper=0.81, s=0.84, s.lower=0.64, s.upper=0.94, 
                 d=2, n=100, CI=TRUE)
estimateN(count=3, p=res.p["p"], p.lower=res.p["2.5\%"], p.upper=res.p["97.5\%"])

}
\keyword{methods}
\keyword{misc}