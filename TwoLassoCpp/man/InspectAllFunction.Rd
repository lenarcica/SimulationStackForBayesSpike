\name{InspectAllFunction}
\alias{InspectAllFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ L0 Minimizer for Regression Problems. }
\description{
  L0 Minimizer Calculation.
}
\usage{
InspectAllFunction(YYtake, XXtake, PenKeep = 1, NKeep = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{YYtake}{ yy Insert Dataset }
  \item{XXtake}{ xx Insert Dataset }
  \item{PenKeep}{ Penalty for every active coefficients }
  \item{NKeep}{ Asks maximum number of covariates to keep. }
}
\details{
  Calculate L0 minimizing Beta Estimate
}
\value{
  Returns Fit Object
}
\references{ Alan Lenarcic's 2009 Thesis }
\author{ Alan Lenarcic }
\note{  
}
\seealso{ EMLARS, EMRIDGE, CoordinateDescent, XLOperation }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
