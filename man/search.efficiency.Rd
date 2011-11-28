\name{search.efficiency}
\alias{search.efficiency}
\title{
Estimates detection probability per person and vegetation density using a binomial model
}
\description{
The function estimates detection probability per vegetation density and person. 
A binomial mixed model is used with vegetation density as fixed effect and person as random factor, when data of more than 2 persons are available.
The number of found and the number of overseen items have to be provided per person and vegetation density class.
The detection probabilities are given per person with a 95 percent credible interval. An average vegetation density specific detection probabilty over all persons is given in addition. 
}
\usage{
search.efficiency(dat, nsim = 1000)
}
\arguments{
  \item{dat}{
a data.frame containing the following variables:\cr 
  person     : names of the persons who searched\cr
  visibility : visibility class\cr
  detected   : number of detected items\cr
  notdetected: number of not detected items\cr
if no visibility classes are available, the variable visibility should contain only one value, e.g. all rows = 0
}
  \item{nsim}{
the number of simulations from the posterior distributions of the model parameters used to construct the 95 percent credible intervals
}
}
\details{
The function uses the function glmer of the package lme4 and the function sim of the package arm.
}
\value{
A list with two elements of class data.frame
  \item{f.perperson}{a data.frame with the estimated detection probabilities per person and visibility class and its lower and uper limit of the 95 percent credible interval}
  \item{f.average}{a data.frame with the estimated detection probabilities per visibility class averaged over the persons}
}
\references{
Gelman A, Hill J (2007) Data Analysis Using Regression and Multilevel and Hierarchical Models. Cambridge: Cambridge Universtiy Press.

Niermann I, Brinkmann R, Korner-Nievergelt F, Behr O (2011) Systematische Schlagopfersuche - Methodische Rahmenbedingungen, statistische Analyseverfahren und Ergebnisse. In: Brinkmann R, Niermann I, Behr O, editors. Entwicklung von Methoden zur Untersuchung und Reduktion des Kollisionsrisikos von Fledermaeusen an Onshore-Windenergieanlagen Goettingen: Cuvillier Verlag. pp. 40-115.
}
\author{
Fraenzi Korner-Nievergelt
}

\examples{
data(searches)
searches
search.efficiency(searches)

}
\keyword{misc}
