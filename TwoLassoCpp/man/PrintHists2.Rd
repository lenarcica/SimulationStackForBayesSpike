\name{PrintHists2}
\alias{PrintHists2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Old histogram printing for simulations }
\description{
  An Old function no longer in use for printing simulations.
}
\usage{
PrintHists2(GFMA1, TCX = c(3, 4, 6), TCY = c(1, 2, 3), SixTeller = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GFMA1}{ Output from SimForMe2 }
  \item{TCX}{ On X axis, list of estimators from 2-11 desired to compare }
  \item{TCY}{ 1,2,3 refer to various metrics (Type II, Type I, and Beta Differential }
  \item{SixTeller}{ A seemingly unused Flag }
}
\details{
   I believe this is un used anymore.  
}
\value{
  Outputs Nothing
}
\references{ "TwoLasso for Bayes Adapted Selection", Lenarcic, Journal of
   Computational and Graphical Statistics, Submitted 2009. }
\author{ Alan Lenarcic }
\note{  Nothing to be Noted
}
\seealso{ PrintHistsV2, PrintMyPlotter }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ iplot }% __ONLY ONE__ keyword per line
