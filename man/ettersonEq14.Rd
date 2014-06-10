\name{ettersonEq14}
\alias{ettersonEq14}
\title{
Equation 14 of Etterson (2013) Ecological Applications 23, 1915-1925
}
\description{
Calculates the probability that a carcass that has fallen in the search area during 
n regular or irregular search intervals is found  by a searcher. 
}
\usage{
ettersonEq14(s, f, J)
}
\arguments{
  \item{s}{persistence probability, probability that a carcass persists a specific time interval 
  (normally a day), i.e. daily persistence probability. Note that the same time unit 
  as for the length of the search interval must be used.
}
  \item{f}{searcher efficiency, i.e. the probability that a carcass present in the search area
  is found during one search.
}
  \item{J}{a vector containing the lengths of the search intervals.
}
}
\details{
 This formula can handle irregular search intervals. It assumes constant persistence time
 and constant searcher efficiency.  
 Confidence intervals are not yet implemented, but a Monte Carlo confidence intervals 
 will be calculated by the function estimateN with the argument form="etterson". 
}
\value{
the probability that a carcass that has fallen into the search area during the study
is found by the searcher.
}
\references{
Etterson, M.A. (2013) Hidden Markov models for estimating animal mortality from 
antropogenic hazards. Ecological Applications, 23, 1915-1925.
}
\author{
M. A. Etterson
}
\note{
 In case of regular search intervals the functions ettersonEq14 and pkorner 
 are equivalent (see first example). A mathematical proof for this equivalence is
 given in the tutorial to carcass (Korner-Nievergelt et al., re-submitted to Wildlife Biology June 2014).
 When search intervals are irregular, 
 pkorner overestimates the detecton probability (i.e. underestimates fatality rates)
 as you see in the second example below.
}

\seealso{
\code{\link{ettersonEq14v1}}
\code{\link{ettersonEq14v2}}
\code{\link{pkorner}}
\code{\link{phuso}}
\code{\link{perickson}}
}

\examples{

# in case of regular search intervals, the calculations below give the same results
ettersonEq14(s=0.8, f=0.8, J=c(3,3,3,3,3))
pkorner(s=0.8, f=0.8, d=3, n=5)


 # in case of irregular search intervals the function ettersonEq14 is more appropriate
ettersonEq14(s=0.8, f=0.8, J=c(3,5,1,4,2))
pkorner(s=0.8, f=0.8, d=mean(c(3,5,1,4,2)), n=5)

}

\keyword{methods}
\keyword{misc}