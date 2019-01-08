\name{NamingClenture2}
\alias{NamingClenture2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ NamingClenture2 }
\description{
  Names Files from Simulations
}
\usage{
NamingClenture2(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps, BRealInput = c(1, -1, 1, -1, 1, -1), LARSSeekFlag = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{ppiuseV}{ Unused Active variable}
  \item{sigmaFactCoV}{ Mean and variance of active correlations between covariates }
  \item{ppCovV}{ Probability of active correlation between two covariates }
  \item{NNOnV}{ Number of Samples per simulation }
  \item{NROnV}{ Total Number of Covariates }
  \item{NoiseOnV}{ Noise level of data }
  \item{NumReps}{ Number of repetitions (200 default) }
  \item{BRealInput}{ Values of True Active coefficients for beginning }
  \item{LARSSeekFlag}{LARSSeekFlag input value for GetForMeAAA2 }
}
\details{
  Internal Naming function
}
\value{
  A String
}
\references{  }
\author{ Alan Lenarcic }
\note{ Internal Naming Function
}
\seealso{  LoadTwoCC }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps, 
    BRealInput = c(1, -1, 1, -1, 1, -1), LARSSeekFlag = 0) 
{
    if (!is.vector(BRealInput)) {
        UseStart <- "NonV"
    }
    else if (length(BRealInput) == 6 && all(BRealInput == c(1, 
        -1, 1, -1, 1, -1))) {
        UseStart <- InitiateS2
    }
    else if (length(BRealInput) == 4 && all(BRealInput == c(4, 
        3, -2.5, 1))) {
        UseStart <- InitiateS
    }
    else {
        UseStart <- "Alt"
    }
    if (LARSSeekFlag == 0) {
        NC <- paste(UseStart, "A", tSeq(ppiuse), "B", tSeq(sigmaFactCov), 
            "C", tSeq(ppCov), "D", tSeq(NNOn), "E", tSeq(NROn), 
            "F", tSeq(NoiseOn), "G", tSeq(NumReps), ".csv", sep = "")
    }
    else {
        NC <- paste(UseStart, "LS", tSeq(LARSSeekFlag), "A", 
            tSeq(ppiuse), "B", tSeq(sigmaFactCov), "C", tSeq(ppCov), 
            "D", tSeq(NNOn), "E", tSeq(NROn), "F", tSeq(NoiseOn), 
            "G", tSeq(NumReps), ".csv", sep = "")
    }
    return(NC)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ print }% __ONLY ONE__ keyword per line
