\name{GetForMeFunction2}
\alias{GetForMeFunction2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ GetForMeFunction2 }
\description{
    GetForMeFunction2
}
\usage{
GetForMeFunction2(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps, FilesInDIR = -999, PrintFlag = 0, PrintOutFlags = 1, BRealInput = c(1, -1, 1, -1, 1, -1), LARSSeekFlag = 0, musize = 10)
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
  \item{PrintFlag}{ Whether the algorithms should print information during run}
  \item{PrintOutFlags}{ Whether Flags should be printed during run }
  \item{BRealInput}{ The active coefficients of the true Beta vector }
  \item{LARSSeekFlag}{ Whether LARSSeekFlag is a non zero amount }
  \item{musize}{ musize targets prior information for EM2PiLasso }
}
\details{
  Not Much listed
}
\value{

}
\references{ ~put references to the literature/web site here ~ }
\author{ Alan Lenarcic }
\note{
}
\seealso{ GetForMeAAA2, GetForMeAll2 }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ datasets }% __ONLY ONE__ keyword per line
