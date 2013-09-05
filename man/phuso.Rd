\name{phuso}
\alias{phuso}
\title{
Carcass detection probability acording to Huso 2010
}
\description{
Calculates carcass detection probability based on persistence probability, searcher efficiency, search interval and duration of study.
}
\usage{
phuso(s, t.bar, f, d)
}
\arguments{
  \item{s}{
persistence probability: probability that a carcass remains on the study plot for one day. This parameter is only used when t.bar is not provided.
}
  \item{t.bar}{
mean persistence time (in days). Alternatively, s can be provided.  
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
