\name{perickson}
\alias{perickson}
\title{
Carcass detection probability acording to Erickson et al. 2004
}
\description{
Calculates carcass detection probability based on persistence time, searcher efficiency, search interval and duration of study.
}
\usage{
perickson(t.bar=NULL, s, f, d)
}
\arguments{
  \item{t.bar}{
mean carcass persistence time in days (or another time unit)
}
  \item{s}{
mean (e.g. daily) carcass persistence probability. 
}

  \item{f}{
searcher efficiency: probability that a carcass present on the study plot is detected by a searcher during a seach.
}
  \item{d}{
search interval: time (in days) between two searches
}
}
\details{
Time measurements (search interval and carcass persistence time or persistence probability) should be given in the same unit. 
Either t.bar or s has to be provided.
}
\value{
a proportion: the probability that an animal that dies during the study period on the study plot is detected by a searcher
}
\references{
Erickson WP, Jeffrey J, Kronner K, Bay K (2004) Stateline Wind Project Wildlife. Monitoring Final Report, July 2001 - December 2003. - Technical report peer-reviewed by and submitted to FPL Energy, the Oregon Energy Facility Siting Council, and the Stateline Technical Advisory Committee, Oregon, USA, 105 pp
}
\author{
Fraenzi Korner-Nievergelt
}

\seealso{
\code{\link{phuso}}
\code{\link{pkorner}}
}
\examples{
perickson(t.bar=30, f=0.8, d=1)
}
\keyword{methods}
\keyword{misc}
