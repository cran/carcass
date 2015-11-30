\name{search.efficiency}
\alias{search.efficiency}
\title{
Estimates detection probability per person and visibility classes using a binomial model
}
\description{
The function estimates detection probability per visibility class and person.
A binomial mixed model is used with vegetation density as fixed effect and person as random factor, when data of more than 2 persons are available.
The number of found and the number of overseen items have to be provided per person and visibility class (one line per observer and visibility class). When only one visibility class is available, the variable visibility must be provided with only one entry. It does not matter what.
The detection probabilities are given per person with a 95\% credible interval. An average vegetation density specific detection probabilty over all persons is given in addition. 
}
\usage{
search.efficiency(dat=NA, person=NA, visibility=NA, detected=NA,
    notdetected=NA, nsim = 1000)
}
\arguments{
  The search efficiancy data may be provided as a \code{data.frame} containing all data or, alternatively, as seperate vectors. If no visibility classes are available, the variable visibility should be a vector with the same length as the others containing a single value (e.g. "not_measured").
  \item{dat}{\code{Data.frame} containing the following columns:\cr 
  \code{person}: names of the persons who searched\cr
  \code{visibility}: visibility class\cr
  \code{detected}: number of detected items\cr
  \code{notdetected}: number of not detected items\cr}  
  \item{person}{vector with names of the persons who searched}
  \item{visibility}{vector with visibility classes}
  \item{detected}{numeric vector with number of detected items}
  \item{notdetected}{numeric vector with number of not detected items}
  \item{nsim}{
the number of simulations from the posterior distributions of the model parameters used to construct the 95 percent credible intervals
}
}
\details{
The function uses the function glmer of the package lme4 and the function sim of the package arm.
}
\value{
A list with two elements of class data.frame
  \item{f.perperson}{a \code{data.frame} with the estimated detection probabilities per person and visibility class, its lower and uper limit of the 95\% credible interval and its standard error}
  \item{f.average}{a data.frame with the estimated detection probabilities per visibility class averaged over the persons}
}
\references{
Gelman A, Hill J (2007) Data Analysis Using Regression and Multilevel and Hierarchical Models. Cambridge: Cambridge University Press.

Niermann I, Brinkmann R, Korner-Nievergelt F, Behr O (2011) Systematische Schlagopfersuche - Methodische Rahmenbedingungen, statistische Analyseverfahren und Ergebnisse. In: Brinkmann R, Niermann I, Behr O, editors. Entwicklung von Methoden zur Untersuchung und Reduktion des Kollisionsrisikos von Fledermaeusen an Onshore-Windenergieanlagen Goettingen: Cuvillier Verlag. pp. 40-115.
}
\author{
Fraenzi Korner-Nievergelt
}

\examples{
data(searches)
searches

# Call to the function with data provided as data.frame:
\dontrun{search.efficiency(searches)}

# Alternative:
per <- searches$person
visi <- searches$visibility
det <- searches$detected
notdet <- searches$notdetected
\dontrun{search.efficiency(person=per, visibility=visi, detected=det, notdetected=notdet)}

}
\keyword{misc}
