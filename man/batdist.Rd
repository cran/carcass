\name{batdist}
\alias{batdist}
\docType{data}
\title{
Distribution of bat carcasses below wind turbines
}
\description{
The data contains the proportion of bats killed at wind turbines falling into different 10m-distance rings from the turbine.
It is a collection of empirical and theoretical estimates from different publications.  
}
\usage{data(batdist)}
\format{
  A data frame with 6 observations on the following 16 variables.
  \describe{
    \item{\code{studyname}}{a factor with the name of the study}
    \item{\code{study}}{a factor with the code of the study, see tutorial}
    \item{\code{rotor.diameter.m}}{a numeric vector giving the average rotor diameter (in m) for the turbines investigated}
    \item{\code{nacelle.height.m}}{a numeric vector giving the average nacelle height (in m) for the turbines investigated}
    \item{\code{n}}{a numeric vector indicating the sample size (i.e. numbers of carcasses found)}
    \item{\code{method}}{a factor indicating the method used to obtain the carcass distribution. \code{empirical}: based on field data, \code{physical_model}: according to a theoretical ballistic model}
    \item{\code{species}}{a factor with only one level at the moment \code{bat}}
    \item{\code{prop10m}}{a numeric vector with the proportion of killed bats lying in between 0 and 10 m from the turbine}
    \item{\code{prop20m}}{a numeric vector with the proportion of killed bats lying in between 10 and 20 m from the turbine}
    \item{\code{prop30m}}{a numeric vector with the proportion of killed bats lying in between 20 and 30 m from the turbine}
    \item{\code{prop40m}}{a numeric vector with the proportion of killed bats lying in between 30 and 40 m from the turbine}
    \item{\code{prop50m}}{a numeric vector with the proportion of killed bats lying in between 40 and 50 m from the turbine}
    \item{\code{prop60m}}{a numeric vector with the proportion of killed bats lying in between 50 and 60 m from the turbine}
    \item{\code{prop70m}}{a numeric vector with the proportion of killed bats lying in between 60 and 70 m from the turbine}
    \item{\code{prop80m}}{a numeric vector with the proportion of killed bats lying in between 70 and 80 m from the turbine}
    \item{\code{prop50mplus}}{a numeric vector: the sum of prop60m, prop70m and prop80m}
  }
}

\references{
Arnett, E. B., W. P. Erickson, J. Kerns, and J. Horn. 2005. Relationships between Bats and Wind Turbines in Pennsylvania and West Virginia:
An Assessment of Fatality Search Protocols, Patterns of Fatality, and Behavioral Interactions with Wind Turbines. Bat Conservation International, Austin, Texas, USA.

Hull, C. L. and S. Muir. 2010. Search area for monitoring bird and bat carcasses at wind farms using a Monte-Carlo model. Australian Journal of Environmental Management 17:77-87.

Niermann, I., R. Brinkmann, F. Korner-Nievergelt, and O. Behr. 2011. Systematische Schlagopfersuche - Methodische Rahmenbedingungen, statistische Analyseverfahren und Ergebnisse. Pages 40-115 in R. Brinkmann, O. Behr, I. Niermann, and M. Reich, editors. Entwicklung von Methoden zur Untersuchung und Reduktion des Kollisionsrisikos von Fledermaeusen an Onshore-Windenergieanlagen. Cuvillier Verlag, Goettingen.

}
\examples{
data(batdist)
batdist
}
\keyword{datasets}
