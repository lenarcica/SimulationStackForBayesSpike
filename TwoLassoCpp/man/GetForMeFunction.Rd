\name{GetForMeFunction}
\alias{GetForMeFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function for Simulating Data, Old}
\description{
  Old Version of Simulating Data and Estimators
}
\usage{
GetForMeFunction(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps, FilesInDIR = -999, musize = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppiuse}{ An Integer that doesn't have too much a role anymore. }
  \item{sigmaFactCov}{ the variance and mean of non-zero correlations }
  \item{ppCov}{ probability of a nonzero connection appearing between variables }
  \item{NNOn}{ Number of sample size for simulations. }
  \item{NROn}{ Number of variables for simulations }
  \item{NoiseOn}{ Sigma squared value }
  \item{NumReps}{ Number of repetitions per simulation parameters }
  \item{FilesInDIR}{ A file name to save to }
  \item{musize}{ musize targets prior information for EM2PiLasso }
}
\details{

}
\value{

}
\references{Alan Lenarcic 2009 Thesis }
\author{ Alan Lenarcic }
\note{  Not Much
}
\seealso{ GetforMeFunction2 }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ methods }% __ONLY ONE__ keyword per line
