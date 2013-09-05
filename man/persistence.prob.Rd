\name{persistence.prob}
\alias{persistence.prob}
\title{
Estimates carcass persistence probability based on carcass removal experiment data
}
\description{
This function either uses a Cox proportional hazard model or an exponential model (if persistence probability can assumed to be constant over time) to estimate daily persistence probabilities of carcasses. 
}
\usage{
persistence.prob(turbineID, perstime, status, pers.const = FALSE, R = 10000)
}
\arguments{
  \item{turbineID}{
factor, character or numeric vector with name of the turbines or sites for which the carcass persistence probability should be estimated separatedly\cr
}
  \item{perstime}{  
numeric vector with the persistence times for each object (e.g. in days)
}
  \item{status}{
indicator variable of observed removal (1= removal has been observed, 0 = object was still there at the end of the observation period)\cr
}
  \item{pers.const}{
logical value indicating whether a constant persistenc probability over time can be assumed. default is FALSE. If FALSE, a Cox proportional hazard model is used and for each turbine/site the estimated proportion of remaining objects is given for each day. If TRUE, an exponential model is fitted and the estimated daily persistence probability is give for each turbine/site.
}
  \item{R}{
number of Monte Carlo simulations used to obtain the 95 percent confidence intervals of the estimated persistence probabilities from the exponential model
}
}
\details{
Note that there is increasing evidence in the literature that carcass persistence probability increases with the age of a carcass. Thus you are saver to use non-constant persistence probabilities unless you have tested, how seriously an assumption of constant persistence probability influences your results. 
}
\value{
If you do not assume constant persistence probability, the function returns a list:
  \item{persistence.prob}{matrix with estimated proportion of remaining carcasses after each time indicated on the rows for each turbine/site indicated on the column}
  \item{estpers.lwr}{the lower limits of the 95 percent confidence intervals}
  \item{estpers.upr}{the upper limits of the 95 percent confidence intervals}
If you assume constant persistence probability, the function returns a data frame with the following variables:\cr
  turbineID: name of the turbine/site\cr
  persistence.prob: estimated daily persistence probability\cr
  lower: lower limit of the 95 percent confidence interval of the estimated persistence probability\cr
  upper: upper limit of the 95 percent confidence interval of the estimated persistence probability\cr
  mean.persistence.time: estimated mean persistence time\cr 
}
\references{
Cox, D. R. 1972. Regression models and life-tables (with discussion). Journal of the Royal Statistical Society B 34:187-220.\cr
Klein, J. P. and M. L. Moeschberber. 2003. Survival Analysis, Techiques for Censored and Truncated Data. Springer, New York.
}
\author{
Fraenzi Korner-Nievergelt
}
\note{
Whether the models used in this function fits to your data is not in the responsibility of the author!
}

\seealso{
\code{\link{coxph}}
\code{\link{survreg}}
}
\examples{
data(persistence)
persistence.prob(persistence$turbineID, persistence$perstime, persistence$status)
persistence.prob(persistence$turbineID, persistence$perstime, persistence$status, 
    pers.const=TRUE)
}

\keyword{methods}
\keyword{misc}