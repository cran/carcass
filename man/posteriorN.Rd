\name{posteriorN}
\alias{posteriorN}
\alias{posterior.N}
\title{
Posterior distribution of the number of fatalities based on the number of carcasses found and the probability of finding a carcass.
}
\description{
This function gives the posterior distribution of the number of fatalities based on the number of carcasses counted during carcass searches and the probability of detecting a carcass. 
The function uses the theorem of Bayes. 
}
\usage{
posteriorN(p, nf = 0, maxN = 1000, ci.int = 0.95, plot = TRUE,
              dist = FALSE)
}
\arguments{
  \item{p}{
probability of finding a carcass. This probability can be obtained by e.g. pkorner, phuso or perickson.
}
  \item{nf}{
the number of carcasses found
}
  \item{maxN}{
the maximal number of fatalities for which a posterior densitiy should be given
}
  \item{ci.int}{
length of the credible interval that should be given. Default is 0.95.
}
  \item{plot}{
if TRUE the posterior distribution is plotted
}
  \item{dist}{
if TRUE the posterior density is given for all elements in 0:maxN
}
}
\value{
\item{interval}{vector of length 2 containing the lower and upper limit of the interval specified in the argument ci.int}
\item{expected}{median of the posterior distribution of the number of fatalities}
\item{HT.estimate}{Horvitz-Thompson estimate, i.e. nf/p}
\item{pN}{posterior density for each element in 0:maxN number of fatalities}
}
\references{
Korner-Nievergelt F, Korner-Nievergelt P, Behr O, Niermann I, Brinkmann R, Hellriegel B (2011) A new method to determine bird and bat fatality at wind energy turbines from carcass searches. Wildlife Biology 17: 350-363
}
\author{
Fraenzi Korner-Nievergelt
}
\note{
This function cannot handle larger numbers, e.g. maxN > 2000 is not possible.
The function does not take into account the uncertainty of the estimate for the detection probability. If you want do account for this uncertainty, use the function estimateN.
}


\seealso{
\code{\link{estimateN}}
}
\examples{
posteriorN(p=0.5, nf=3, dist=TRUE, maxN=15)
}

\keyword{methods}
\keyword{misc}
