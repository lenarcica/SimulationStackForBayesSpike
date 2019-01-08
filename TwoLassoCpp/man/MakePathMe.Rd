\name{MakePathMe}
\alias{MakePathMe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  Tries to define inital "PathMe" where data is saved.
}
\usage{
MakePathMe()
}
%- maybe also 'usage' for other objects documented here.
\details{
  Used for saving simulations and tables
}
\value{

}
\references{ ~put references to the literature/web site here ~ }
\author{ Alan Lenarcic }
\note{ May have to be fixed for other systems.
}
\seealso{ LoadSavedOutPutDirectory, LoadSavedTableDirectory }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function () 
{
    PathMe <<- paste(.Library, "/TwoLasso/", sep = "")
    return(PathMe)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ documentation }% __ONLY ONE__ keyword per line
