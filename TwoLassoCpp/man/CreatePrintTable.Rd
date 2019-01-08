\name{CreatePrintTable}
\alias{CreatePrintTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print Table Helper }
\description{
    Given Simulation output, summarizes into R table for later print functions.
}
\usage{
CreatePrintTable(GFMAA, OneVV, KPAm, PrMeVec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GFMAA}{ Table of SimForMe2 Output }
  \item{OneVV}{ Matrix with rows of Parameter settings }
  \item{KPAm}{ Number of ActiveCovariates in Set }
  \item{PrMeVec}{ List of subset integers 1-10 for the ten different estimators}
}
\details{
  Plotting Function
}
\value{
  Returns that Table
}
\references{ Lenarcic 2009 Thesis }
\author{ Alan Lenarcic }
\note{  Helper Function
}
\seealso{ DoAllTheSaving  }
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ distribution }% __ONLY ONE__ keyword per line
