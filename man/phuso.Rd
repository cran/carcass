\name{phuso}
\alias{phuso}
\title{
Carcass detection probability acording to Huso 2010
}
\description{
Calculates carcass detection probability based on persistence time, searcher 
efficiency, search interval and duration of study. Persistence time 
obtained by any survival model (e.g. Exponential, Weibull, Log-Normal,...) can
be given (but see below). Alternatively to persistence time, persistence probability can be given,
In this case, the exponential model (constant persistence probability) is used to transform
persistence probability into mean persistence time. \cr
An argument shape and distribution for other than the Exponential persistence model will be provided soon.
}
\usage{
phuso(s, t.bar, f, d)
}
\arguments{
  \item{s}{
persistence probability: probability that a carcass remains on the study plot for one day. 
This parameter is only used when t.bar is not provided. Note that if s is given instead
of t.bar, an exponential survival function (constant persistence probability) is used.
}
  \item{t.bar}{
mean persistence time (in days). Alternatively, s can be provided. Then, an exponential persistence model is assumed.  
}
  \item{f}{
searcher efficiency: probability that a carcass present on the study plot is detected by a searcher during a seach.
}
  \item{d}{
search interval: time (in days) between two searches
}
}
\details{
Time measurements (search interval) and reference time units (persistence probability) should be given in the same unit. 
}
\value{
a proportion: the probability that an animal that dies during the study period on the study plot is detected by a searcher
}
\references{
Huso M (2010) An estimator of wildlife fatality from observed carcasses. Environmetrics 22: 318-329
}
\author{
Fraenzi Korner-Nievergelt
}

\seealso{
\code{\link{perickson}}
\code{\link{pkorner}}
}
\examples{
phuso(s=0.8, f=0.7, d=7)
}
\keyword{methods}
\keyword{misc}
