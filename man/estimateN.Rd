\name{estimateN}
\alias{estimateN}
\title{Estimation of number of killed animals based on carcass searches and estimates for detection probability
}
\description{
Estimates the number of killed animals from carcass searches taking into account the uncertainty in the estimates for carcass detection probability (or alternatively for searcher efficiency and carcass persistence probability).
}
\usage{
estimateN(count, p = NA, p.lower = NA, p.upper = NA,
          f = NA, f.lower = NA, f.upper = NA, s = NA, 
          s.lower = NA, s.upper = NA, pform = "korner", d = 1, 
          n = NA, maxn = 500, nsim = 1000, plot = TRUE, postdist=FALSE)
}
\arguments{
  \item{count}{
number of carcasses found
}
  \item{p}{
estimate for detection probability (alternatively give f, s, d and n)
}
  \item{p.lower}{
lower limit of 95percent CI of detection probability
}
  \item{p.upper}{
upper limit of 95percent CI of detection probability
}
  \item{f}{
estimate for searcher efficiency
}
  \item{f.lower}{
lower limit of 95percent CI of searcher efficiency
}
  \item{f.upper}{
upper limit of 95percent CI of searcher efficiency
}
  \item{s}{
estimate for persistence probability
}
  \item{s.lower}{
lower limit of 95percent CI of persistence probability
}
  \item{s.upper}{
upper limit of 95percent CI of persistence probability
}
  \item{pform}{
formula used to estimate p, one of "korner", "huso", "erickson"
}
  \item{d}{
 search interval, the number days (or any other time unit) between two searches
}
  \item{n}{
number of searches
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
}
\details{
A documented code of this function if given in the Appendix of Korner-Nievergelt et al. (2011) J. Wildlife. Biol.

If postdist is TRUE, the output is a list, if postdist is FALSE, the output is a named numeric vector. 
}
\value{
Either a list or a numeric vector containing the estimate (median of posterior) and the 95percent credible interval
}
\references{
Erickson WP, Jeffrey J, Kronner K, Bay K (2004) Stateline Wind Project Wildlife. Monitoring Final Report, July 2001 - December 2003. - Technical report peer-reviewed by and submitted to FPL Energy, the Oregon Energy Facility Siting Council, and the Stateline Technical Advisory Committee, Oregon, USA, 105 pp\cr
Huso M (2010) An estimator of wildlife fatality from observed carcasses. Environmetrics 22: 318-329\cr
Korner-Nievergelt F, Korner-Nievergelt P, Behr O, Niermann I, Brinkmann R, Hellriegel B (2011) A new method to determine bird and bat fatality at wind energy turbines from carcass searches. Wildlife Biology 17: 350-363
}
\author{
F. Korner-Nievergelt
}


\seealso{
\code{\link{posterior.N}}
}
\examples{

estimateN(count=3, f=0.72, f.lower=0.62, f.upper=0.81, s=0.84, s.lower=0.64, 
            s.upper=0.94, d=2, pform="korner", n=100, maxn=500, nsim=1000, 
            plot=TRUE)

estimateN(count=3, f=0.72, f.lower=0.62, f.upper=0.81, s=0.84, s.lower=0.64, 
            s.upper=0.94, d=2, pform="huso", maxn=500, nsim=1000, plot=TRUE)

estimateN(count=3, p=0.6, p.lower=0.5, p.upper=0.7)
}
\keyword{methods}
\keyword{misc}