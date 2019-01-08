\name{GiveGiveMeDiffD}
\alias{GiveGiveMeDiffD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ GiveGiveMe Difference Distances }
\description{
  Used for comparing difference in magnitude of setimates
}
\usage{
GiveGiveMeDiffD(ArrsCoefs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ArrsCoefs}{ Matrix of Beta Estimates for different estimators. }
}
\details{
  Calculates sum squared deviance between variables
}
\value{
  Returns vector long list of distances
}
\references{ Alan Lenarcic's Harvard Thesis 2009 }
\author{ Alan Lenarcic }
\note{ 
}
\seealso{ GiveGiveMe }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ArrsCoefs) 
{
    ALT <- length(ArrsCoefs[, 1])
    RetAD <- vector("numeric", ALT * (ALT - 1)/2)
    ndn <- 1
    for (ii in 1:(ALT - 1)) {
        for (jj in (ii + 1):ALT) {
            AD1 <- ArrsCoefs[ii, ]
            AD2 <- ArrsCoefs[jj, ]
            ADD <- AD1 - AD2
            if ((is.nan(AD1[1])) || is.nan(AD2[1]) || AD1[1] < 
                -999 || AD2[1] < -999) {
                RetAD[ndn] <- -10
            }
            else {
                RetAD[ndn] <- sum(ADD^2)
            }
            ndn <- ndn + 1
        }
    }
    return(RetAD)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ documentation }% __ONLY ONE__ keyword per line
