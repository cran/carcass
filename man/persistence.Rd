\name{persistence}
\alias{persistence}
\docType{data}
\title{
Times until removal for brown mice carcasses 
}
\description{
Experimental data to measure persistence time of bat-sized mammal carcasses. In total, 595 brown mice and 41 bat carcasses have been
put on the ground at 30 different sites (turbineID). During 14 days, the carcasses have been observed and the number of days until their
disappearance was noted. All carcasses that did not disappear until day 14 were marked as censored (status=0) in the data.  
}
\usage{data(persistence)}
\format{
  A data frame with 636 observations on the following 3 variables.
  \describe{
    \item{\code{turbineID}}{a factor with site or turbine names}
    \item{\code{perstime}}{a numeric vector containing the days until disappearance (persistence times)}
    \item{\code{status}}{a numeric vector 1 = disappearance is observed, 0 = no disappearance until end of study (censored)}
  }
}
\details{
Due to conservation reasons, experiments to estimate carcass persistence times have been conducted with brown mice instead of bat carcasses. 
There has been concern that bats have longer persistence times than mice. In our data, during the first four days removal probability of mice and bats were identical. 
After four days, bats seem to be removed less likely. However, this difference was based on such low sample size that we do not dare to claim any difference at the moment.
For the moment, we can say that it is probably save to use mice when the search interval is less than four days. But for longer search intervals, it may be better to use bats to measure carcass persistence probability.\cr
}
\source{
Niermann, I., R. Brinkmann, F. Korner-Nievergelt, and O. Behr. 2011. Systematische Schlagopfersuche - Methodische Rahmenbedingungen, statistische Analyseverfahren und Ergebnisse. Pages 40-115 in R. Brinkmann, O. Behr, I. Niermann, and M. Reich, editors. Entwicklung von Methoden zur Untersuchung und Reduktion des Kollisionsrisikos von Fledermaeusen an Onshore-Windenergieanlagen. Cuvillier Verlag, Goettingen.
}
\references{
Niermann, I., R. Brinkmann, F. Korner-Nievergelt, and O. Behr. 2011. Systematische Schlagopfersuche - Methodische Rahmenbedingungen, statistische Analyseverfahren und Ergebnisse. Pages 40-115 in R. Brinkmann, O. Behr, I. Niermann, and M. Reich, editors. Entwicklung von Methoden zur Untersuchung und Reduktion des Kollisionsrisikos von Fledermaeusen an Onshore-Windenergieanlagen. Cuvillier Verlag, Goettingen.

}
\examples{
data(persistence)
head(persistence)
}
\keyword{datasets}
