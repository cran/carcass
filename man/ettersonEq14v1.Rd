\name{ettersonEq14v1}
\alias{ettersonEq14v1}
\title{
Equation 14 of Etterson (2013) Ecological Applications 23, 1915-1925, 
adapted so that persistence probability and searcher efficiency can vary 
with calender date
}
\description{
Calculates the probability that a carcass that has fallen in the search area during 
n regular or irregular search intervals is found  by a searcher.
}
\usage{
ettersonEq14v1(s, f, J)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{vector of persistence probabilities with N elements, where N is the number of days in
the full study period (i.e, N = sum(J)).
}
  \item{f}{vector of searcher efficiency values (probabilities that a carcass 
  present in the search area is found during one search) with n elements, 
  where n is the number of searches performed.
}
  \item{J}{a vector containing the lengths of the search intervals.
}
}
\details{
This function assumes that persistence probability and searcher efficiency change over
time and that this change depends upon serial date within your study period. In other words, it
assumes that all carcasses in the environment on a given date experience the same scavenging probability
and detectability, regardless of how old the carcasses is.
However it does allow persistence probability and searcher efficiency to change with the calendar date.
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
\code{\link{ettersonEq14v2}}
\code{\link{pkorner}}
\code{\link{phuso}}
\code{\link{perickson}}
}
\examples{

J <- c(2,3,2,4,3,5,3,2,3,4)
s <- plogis(seq(0.2, 2, length=sum(J)))
f <- plogis(seq(1.5, 0.9, length=length(J)))
ettersonEq14v1(s,f,J)

}

\keyword{methods}
\keyword{misc}