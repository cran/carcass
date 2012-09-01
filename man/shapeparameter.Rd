\name{shapeparameter}
\alias{shapeparameter}
\title{
Shapeparameters of a beta-distribution from the mean, the lower and upper limit of the 95percent confidence or credible interval
}
\description{
The function returens the shapeparameters a ad b of a beta-distribution from the mean, the lower and upper limit of the 95percent confidence or credible interval
}
\usage{
shapeparameter(m, lwr, upr)
}

\arguments{
  \item{m}{
mean
}
  \item{lwr}{
lower limit of the 95 percent confidence/credible interval
}
  \item{upr}{
upper limit of the 95 percent confidence/credible interval
}
}
\value{
\item{a}{parameter a of a beta-distribution}
\item{b}{parameter b of a beta-distribution}
}
\references{
no reference, i.e. developed from the definition of a and b in a beta distribution obtained e.g. from Gelman et al. (2004)
}
\author{
Fraenzi Korner-Nievergelt
}

\examples{

a <- shapeparameter(0.8, 0.72, 0.88)$a
b <- shapeparameter(0.8, 0.72, 0.88)$b
x <- seq(0, 1, by=0.01)
y <- dbeta(x, a, b)
plot(x, y, type="l")

}
\keyword{methods}
\keyword{misc}
