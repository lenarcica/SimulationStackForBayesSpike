\name{NamingClenture}
\alias{NamingClenture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Old Function to name simulation files }
\description{
  Names Files that are saved for database purposes
}
\usage{
NamingClenture(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps)
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
}
\details{
  Returns a character String for Printing
}
\value{
 A string
}
\references{ }
\author{ Alan Lenarcic }
\note{  Not Really Used for Naming Files anymore
}
\seealso{ NamingClenture2 }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps) 
{
    NC <- paste(InitiateS, "A", tSeq(ppiuse), "B", tSeq(sigmaFactCov), 
        "C", tSeq(ppCov), "D", tSeq(NNOn), "E", tSeq(NROn), "F", 
        tSeq(NoiseOn), "G", tSeq(NumReps), ".csv", sep = "")
    return(NC)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ documentation }% __ONLY ONE__ keyword per line
