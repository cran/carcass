\name{pkorner}
\alias{pkorner}                                                 
\alias{pkorner.version2}                                                 
\title{
Carcass detection probability acording to Korner-Nievergelt et al. 2011
}
\description{
Calculates carcass detection probability based on persistence probability, searcher efficiency, search interval and duration of study.
}
\usage{
pkorner(s, f, d, n, k=0.25, search.efficiency.constant=TRUE)
}
\arguments{
  \item{s}{
persistence probability: probability that a carcass remains on the study plot for one day (or another time unit).
}
  \item{f}{
searcher efficiency: probability that a carcass present on the study plot is detected by a searcher during a seach.
}
  \item{d}{
search interval: time (in days) between two searches
}
  \item{n}{
number of searches: n*d = duration of study
}
  \item{k}{
factor by which the search efficiency is multiplied after each search for a specific carcass (see details)
}
  \item{search.efficiency.constant}{
TRUE if search efficiency is constant over time, FALSE if it decreases by a factor k with each search.
}


}
\details{
Time measurements (search interval) and reference time units (persistence probability) should be given in the same unit.
The argument k together with search.efficiency.constant=FALSE can be given to account for the fact that during the first search 
carcasses with high detectability are found overproportionatedly, whereas carcasses hard to find remain on the study plot. Thus with each
search, searcher efficiency decreases by the factor k. 
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
\code{\link{phuso}}
\code{\link{perickson}}
}
\examples{

f <- 0.72
s <- 0.84
pkorner(s=s, f=f, d=1, n=100)
pkorner(s=s, f=f, d=7, n=15)
pkorner(s=s, f=f, d=1, n=100, search.efficiency.constant=FALSE)


}
\keyword{methods}
\keyword{misc}
