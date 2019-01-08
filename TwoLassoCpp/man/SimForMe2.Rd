\name{SimForMe2}
\alias{SimForMe2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ SimForMe2 }
\description{
  Simulates multiple fits of one set of parameters to compare estimators
}
\usage{
SimForMe2(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps, musize = 10, SigmaEta = 2, SigmaBar = -999, InputTable = NULL, FileName = NULL, PrintOutFlags = 0, BRealInput = c(1, -1, 1, -1, 1, -1), LARSSeekFlag = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppiuse}{ piA to use for all estimates, usually not configured }
  \item{sigmaFactCov}{ For adding covariance to X matrix }
  \item{ppCov}{ Probability of covariance between X_j and X_j' vector }
  \item{NNOn}{ Sample Size }
  \item{NROn}{ Total number of Covariates }
  \item{NoiseOn}{ ~~Describe \code{NoiseOn} here~~ }
  \item{NumReps}{ ~~Describe \code{NumReps} here~~ }
  \item{musize}{ ~~Describe \code{musize} here~~ }
  \item{SigmaEta}{ ~~Describe \code{SigmaEta} here~~ }
  \item{SigmaBar}{ ~~Describe \code{SigmaBar} here~~ }
  \item{InputTable}{ ~~Describe \code{InputTable} here~~ }
  \item{FileName}{ ~~Describe \code{FileName} here~~ }
  \item{PrintOutFlags}{ ~~Describe \code{PrintOutFlags} here~~ }
  \item{BRealInput}{ ~~Describe \code{BRealInput} here~~ }
  \item{LARSSeekFlag}{ ~~Describe \code{LARSSeekFlag} here~~ }
}
\details{
  Simulate and Fit estimators of many TwoLasso and Related Types
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ Alan Lenarcic }
\note{ Simulates multiple fits of one set of parameters to compare estimators
}
\seealso{ GetForMeAAA2, GetForMeAll2 }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
ppiuse = 7; sigmaFactCov = .2; ppCov = .3; NNOn = 100; NROn= 100; NoiseOn = 1;
NumReps = 10; musize = 10; SigmaEta = 2; SigmaBar = -999; InputTable = NULL;
FileName = NULL; PrintOutFlags = 4;
BRealInput = c(1,-1,1,-1,1,-1);  LARSSeekFlag = 2.5;
SimForMe2(ppiuse = ppiuse, sigmaFactCov=sigmaFactCov, ppCov = ppCov, 
   NNOn = NNOn, NROn = NROn, NoiseOn = NoiseOn, NumReps = NumReps, 
   musize = musize, SigmaEta = SigmaEta, SigmaBar = Sigmabar, 
   InputTable = InputTable, FileName = FileName, 
   PrintOutFlags = PrintOutFlags, 
   BRealInput = c(1, -1, 1, -1, 1, -1), LARSSeekFlag = 2.5);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
