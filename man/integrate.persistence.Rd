\name{integrate.persistence}
\alias{integrate.persistence}

\title{
Integrate persistence probability over the discrete time intervals (e.g. days) 
to account for constant arrival of carcasses 
}
\description{
The function integrates the proportion of carcasses remaining on the ground (measured
at discrete time intervals after death) over the last time interval. This gives for each
discrete time point the proportion of carcasses that have arrived continuously 
(uniformly distributed) during the first time interval. This takes into account
that carcasses of the same cohort (arrived during the same time interval) can 
be of different age (max. age difference is length of time interval, e.g one day in most studies).
 
}
\usage{
integrate.persistence(s, n = 20, d = 1)
}

\arguments{
  \item{s}{
Either a number or a numeric vector. If a number is given (which should be the probability that a 
carcass persists one time interval, e.g. 24 hours), an exponential survivor
function (constant persistence probability) is assumed and the proportion of carcasses
persisting until the end of every time interval is calculated. Thus, the function returns a vector. 
The first number is the proportion of carcasses remaining until the end of the first time interval,
which is the time interval during which the carcasses have continuously arrived at a constant rate. 
The second number gives the proportion of carcasses that have arrived during the first time interval 
and that are still persisting until the end of the second time interval, a.s.o. This proportion is calculated
up to the n-th time interval.
If a vector is given, this vector can be either the proportion of carcasses remaining after 
every time interval after death ("survivor function") or the age dependent probability 
to persist one time interval (e.g. daily persistence probability). In the first case, 
the function returns the proportion of persisting carcasses after each time interval
assuming that they arrived at a constant rate during the first time interval. In 
the second case the same integration leads to a average of the persistence probability
at time t and t-1 and can thus be seen as a correction for the different ages of the
carcasses belonging to the same cohort.
ATTENTION: if you use an exponential survival model (constant persistence probability) together
with the function ettersonEqv2 that requires daily persistence probabilities rather than survivor 
functions, you must not provide s as a single number here. Use rep(s, n) in this case! 
}
  \item{n}{
number of time intervals for which the proportion of persisting carcasses is calculated, default is 20
}
  \item{d}{
This argument is only needed by the function estimateN, where n is the number of searches and
d the search interval. Therefore, the number of days of the study is n*d and it is reasonable
to calcualte the persistence probabilities for a maximum time span that equals the study length. 
}
}
\details{
IMPORTANT: if you set the argument "arrival='uniform'" in the estimateN function, you must 
not use the function integrate.persistence prior to use the estimateN function. Since
estimateN does this integration when "arrival='uniform'". But the integrate.persistence
function is useful in combination with the function pkorner, phuso, perikson, and all
etterson-formulas.
}
\value{
A vector containing the proportion of carcasses of one cohort of carcasses (that
have arrived at a constant rate during the first time interval) that persist until the
end of n time intervals. 
}
\references{
Software in prep by Dan Dalthorp and Manuela Huso
}
\author{
F. Korner
}


\seealso{
\code{\link{estimateN}}
}
\examples{
integrate.persistence(0.5)
integrate.persistence(c(0.8,0.7,0.6,0.55))
}

\keyword{misc}

