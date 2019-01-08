\name{PrintHistsV2}
\alias{PrintHistsV2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print Histograms }
\description{
  Used to print Histograms
}
\usage{
PrintHistsV2(GFMA1, TCX = c(3, 4, 6), TCY = c(1, 2, 3), SixTeller = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{GFMA1}{ Output from SimForMe2 }
  \item{TCX}{ On X axis, list of estimators from 2-11 desired to compare }
  \item{TCY}{ 1,2,3 refer to various metrics (Type II, Type I, and Beta Differential }
  \item{SixTeller}{ A seemingly unused Flag }
}
\details{
  No Details other than function that prints out histograms
}
\value{
  A plot
}
\references{ "TwoLasso for Bayes Adapted Selection", Lenarcic, Journal of
   Computational and Graphical Statistics, Submitted 2009. }
\author{ Alan Lenarcic }
\note{ 
  Not Really
}
\seealso{ PrintHists2, PrintMyPlotter }
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ hplot }% __ONLY ONE__ keyword per line
