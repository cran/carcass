\name{ettersonEq14v2}
\alias{ettersonEq14v2}
\title{
Equation 14 of Etterson (2013) Ecological Applications 23, 1915-1925, 
adapted so that persistence probability and searcher efficiency can vary 
with age of the carcass
}
\description{
Calculates the probability that a carcass that has fallen in the search area during 
n regular or irregular search intervals is found  by a searcher.
}
\usage{
ettersonEq14v2(s, f, J)
}
\arguments{
  \item{s}{vector of persistence probabilities with N elements, where N is the number of days in
the full study period (i.e, N = sum(J)) (so N is the maximum number of days that
a carcass could be present in the environment).
}
  \item{f}{vector of discovery probabilities with n elements, where n is the number of searches
performed.
}
  \item{J}{a vector containing the lengths of the search intervals.
}
}
\details{
This function assumes that persistence probability and searcher efficiency change over
time and that this change depends upon the time elapsed since the carcass was killed (carcass age). 
If you assume that the change in detectability of a carcass is due to its age, 
then this function may not be reliable because it indexes the searcher efficiency 
to the number of searches the carcass has been exposed to.
However, if you assume that searcher efficiency changes with the number of
searches a carcass was exposed to because the easy carcasses are found first, 
then this function should perform well.
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

\seealso{
\code{\link{ettersonEq14}}
\code{\link{ettersonEq14v1}}
\code{\link{pkorner}}
\code{\link{phuso}}
\code{\link{perickson}}
}
\examples{
J <- c(2,3,2,4,3,5,3,2,3,4)
s <- plogis(seq(0.2, 2, length=sum(J)))
f <- plogis(seq(1.5, 0.9, length=length(J)))
ettersonEq14v2(s,f,J)
}

\keyword{methods}
\keyword{misc}