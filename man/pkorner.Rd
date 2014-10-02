\name{pkorner}
\alias{pkorner}                                                                                                
\title{
Carcass detection probability acording to Korner-Nievergelt et al. 2011
}
\description{
Calculates carcass detection probability based on persistence probability, searcher efficiency, search interval and duration of study.
}
\usage{
pkorner(s, s.lower=NA, s.upper=NA, f, f.lower=NA, f.upper=NA, d, n, k=0.25, 
        search.efficiency.constant=TRUE, CI=FALSE, nsim=1000)
}
\arguments{
  \item{s}{
persistence probability: Probability that a carcass remains on the study plot for one day (or another time unit). 
Scalar in case of constant persistence probability, vector with probabilities to 
persist until each day after death ('survivor function', i.e. proportion of remaining carcasses for
each day after death) in case of non-constant persistence probabilities.
}
  \item{s.lower}{
lower limit of 95\% CI of persistence probability
}
  \item{s.upper}{
upper limit of 95\% CI of persistence probability
}
  \item{f}{
searcher efficiency: Probability that a carcass present on the study plot is detected by a searcher during a seach.
}
  \item{f.lower}{
lower limit of 95\% CI of searcher efficiency
}
  \item{f.upper}{
upper limit of 95\% CI of searcher efficiency
}
  \item{d}{
search interval: time (in days) between two searches
}
  \item{n}{
number of searches: n*d = duration of study
}
  \item{k}{
Factor by which the search efficiency is multiplied after each search for a specific carcass (see details).
}
  \item{search.efficiency.constant}{
TRUE if search efficiency is constant over time, FALSE if it decreases by a factor k with each search.
}
  \item{CI}{
TRUE if 95\% CI should be calculated based on a beta distribution and Monte Carlo simulations, FALSE otherwise.
}
  \item{nsim}{
number of Monte Carlo simulations
}

}
\details{
Time measurements (search interval) and reference time units (persistence probability) should be given in the same unit.
The argument k together with search.efficiency.constant=FALSE can be given to account for the fact that during the first search 
carcasses with high detectability are found overproportionatedly, whereas carcasses hard to find remain on the study plot. Thus with each
search, searcher efficiency decreases by the factor k. If search.efficiency.constant=TRUE the value given at the argument k is ignored (i.e. k=1). 
If one assumes that persistence probability is not constant, the function \code{\link[carcass]{persistence.prob}} could be used to provide the data for \code{s} (see examples below).
An estimate for s should be provided for every day for a sufficient time span, e.g. 14 days. Missing estimates for intermediate days have to be interpolated. 
After the last day with an s provided, the function assumes that the proportion of persisting carcasses is zero, e.g. s is assumed to be zero from day 15 and onwards if s is a vector of length 14.
}
\value{
a proportion: the probability that an animal that dies during the study period on the study plot is detected by a searcher
}
\references{
Korner-Nievergelt F, Korner-Nievergelt P, Behr O, Niermann I, Brinkmann R, Hellriegel B (2011) A new method to determine bird and bat fatality at wind energy turbines from carcass searches. Wildlife Biology 17: 350-363
}
\author{
Fraenzi Korner-Nievergelt
}

\seealso{
\code{\link{phuso}},
\code{\link{perickson}},
\code{\link{persistence.prob}}
}
\examples{

### Data
  f <- 0.72
  s <- 0.8
  data(persistence)
  attach(persistence)
  sv <- persistence.prob(turbineID, perstime, status)$persistence.prob[,1]
  sv.lower <- persistence.prob(turbineID, perstime, status)$lower[,1]
  sv.upper <- persistence.prob(turbineID, perstime, status)$upper[,1]
  n <- 4
  d <- 3
  
### Constant search efficiency and constant persistence probability
  pkorner(s=s, f=f, d=d, n=n)
  pkorner(s=s, s.lower=0.6, s.upper=0.9, f=f, f.lower=0.6, f.upper=0.8, 
    d=d, n=n, CI=TRUE)

### Decreasing search efficiency and constant persistence probability
  pkorner(s=s, f=f, d=d, n=n, k=0.25, search.efficiency.constant=FALSE)

### Constant search efficiency and decreasing persistence probability
  pkorner(s=sv, f=f, d=d, n=n)

### Decreasing search efficiency and decreasing persistence probability
  pkorner(s=sv, f=f, d=d, n=n, search.efficiency.constant=FALSE)
  pkorner(s=sv, s.lower=sv.lower, s.upper=sv.upper, f=f, f.lower=0.6, 
    f.upper=0.8, d=d, n=n, search.efficiency.constant=FALSE, CI=TRUE)

}
\keyword{methods}
\keyword{misc}
