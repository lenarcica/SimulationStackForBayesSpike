\name{StandardizeXX}
\alias{StandardizeXX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Standardize the coefficients }
\description{
  Standardizes a covariate or response data item
}
\usage{
StandardizeXX(XX)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XX}{ matrix or column to standardize }
}
\details{
  Standardizes covariate set.
}
\value{
  Returns standardized equivalent of input.
}
\references{ Alan Lenarcic's 2009 Harvard Thesis Chapter 2

   "TwoLasso for Bayes Adapted Selection", Lenarcic, Journal of
   Computational and Graphical Statistics, Submitted 2009.  }
\author{ Alan Lenarcic }
\note{  Nothing Interesting to follow
}
\seealso{ standardize }
\examples{
StandardizeXX( t(1:5+ matrix(rnorm(5 * 20,0,1), 5, 20)) );

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (XX) 
{
    if (length(dim(XX)) == 2) {
        XX <- t((t(XX) - meanC(XX))/sd(XX))
    }
    else {
        XX <- (XX - mean(XX))/sd(XX)
    }
    return(XX)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ arith }% __ONLY ONE__ keyword per line
