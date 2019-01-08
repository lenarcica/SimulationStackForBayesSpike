\name{SignZero}
\alias{SignZero}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Rounded Sign function }
\description{
  Rounded Sign function on data vector.
}
\usage{
SignZero(aa, CZ)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aa}{ Number to check for sign }
  \item{CZ}{ minimum value to be considered zero }
}
\details{
   Internal "sign(x)" function, uses fuzzy rounding for near zero values.
}
\value{
 Outputs vector of -1,0,and 1's.
}
\references{ "TwoLasso for Bayes Adapted Selection", Lenarcic, Journal of
   Computational and Graphical Statistics, Submitted 2009.  }
\author{ Alan Lenarcic }
\note{ Not really
}
\seealso{ sign }
\examples{
  SignZero(1.5, .001);
  SignZero(.1, .001);
  SignZero(.00001, .001);
  SignZero(-40, .001);
## The function is currently defined as
function (aa, CZ) 
{
    raa <- aa * 0
    raa[aa > abs(CZ)] <- 1
    raa[aa < -abs(CZ)] <- -1
    return(raa)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ arith }% __ONLY ONE__ keyword per line
