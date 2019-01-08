\name{tSeq}
\alias{tSeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Modifies a number for printing }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tSeq(Number)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Number}{ ~~Describe \code{Number} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{Internal Naming function }
\author{ Alan Lenarcic }
\note{ None Needed
}
\seealso{ NamingClenture2, NamingClenture }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Number) 
{
    if (Number - floor(Number) > 0) {
        NNR <- 10 * (Number - floor(Number))
        if (NNR - floor(NNR) == 0) {
            return(paste(floor(Number), "p", NNR, sep = ""))
        }
        else {
            NNR <- 100 * (Number - floor(Number))
            if (NNR - floor(NNR) == 0) {
                return(paste(floor(Number), "p", NNR, sep = ""))
            }
            else {
                NNR <- floor(1000 * (Number - floor(Number)))
                return(paste(floor(Number), "p", NNR, sep = ""))
            }
        }
    }
    else {
        return(paste(Number, sep = ""))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ IO }% __ONLY ONE__ keyword per line
