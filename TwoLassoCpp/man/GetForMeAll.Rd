\name{GetForMeAll}
\alias{GetForMeAll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ GetForMeAll }
\description{
  Old Version For Simulating Estimator Performance
}
\usage{
GetForMeAll(ppiuseV, sigmaFactV, ppCovV, NNOnV, NROnV, NoiseOnV, NumReps, PrintFlag = 0, musize = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{ppiuseV}{ Unused Active variable}
  \item{sigmaFactV}{ Mean and variance of active correlations between covariates }
  \item{ppCovV}{ Probability of active correlation between two covariates }
  \item{NNOnV}{ Number of Samples per simulation }
  \item{NROnV}{ Total Number of Covariates }
  \item{NoiseOnV}{ Noise level of data }
  \item{NumReps}{ Number of repetitions (200 default) }
  \item{PrintFlag}{ Whether to print information while simulating }
  \item{musize}{ Length of true active vector. }
}
\details{
  Historical Simulation Function
}
\value{
   Historical Simulation Function
}
\references{ TwoLasso For Bayes Adapted Selection }
\author{ Alan Lenarcic }
\note{ 
}
\seealso{ GetForMeAll2, GetForMeAAA, GetForMeAAA2 }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }
\keyword{ datagen }% __ONLY ONE__ keyword per line
\keyword{ internal }
