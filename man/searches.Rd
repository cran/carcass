\name{searches}
\alias{searches}
\docType{data}
\title{
Data of a searcher efficiency trial
}
\description{
Data containing the number of seen and not-seen items of a searcher efficiency trial performed in 3 different vegetation densities (visibilities) and for various persons.
}
\usage{data(searches)}
\format{
  A data frame with 17 observations on the following 4 variables.
  \describe{
    \item{\code{person}}{a factor with the names of the persons}
    \item{\code{visibility}}{a factor with visibility classes (vegetation densities)}
    \item{\code{detected}}{a numeric vector indicating the number of found items}
    \item{\code{notdetected}}{a numeric vector indicating the number of not found items}
  }
}
\source{
Niermann I, Brinkmann R, Korner-Nievergelt F, Behr O (2011) Systematische Schlagopfersuche - Methodische Rahmenbedingungen, statistische Analyseverfahren und Ergebnisse. \emph{In:} Brinkmann R, Niermann I, Behr O, editors. Entwicklung von Methoden zur Untersuchung und Reduktion des Kollisionsrisikos von Fledermaeusen an Onshore-Windenergieanlagen Goettingen: Cuvillier Verlag. pp. 40-115.
}
\examples{
data(searches)
str(searches)
search.efficiency(searches)
}
\keyword{datasets}
