\name{MyLarsGiveMeClose}
\alias{MyLarsGiveMeClose}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ LarsGiveMeClose Again }
\description{
  A program for seeking minimum branching point for "IntClose" active factors
     in a Lars search path.
}
\usage{
MyLarsGiveMeClose(MyLarsOb, IntClose, MinMin = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MyLarsOb}{ Output object from function lars }
  \item{IntClose}{ Number of active factors desired }
  \item{MinMin}{ Minimum Lasso Lambda parameter value. }
}
\details{
}
\value{
  Outputs an integer of the desired step in LARS path with a certain number of
    active factors;
}
\references{ Efron 2004, lars package }
\author{ Alan Lenarcic }
\note{ 
}
\seealso{ lars, LarsCC2, CoordinateDescent }
\examples{
   library(lars); data(diabetes);
	    HitLars2Ob <- lars(diabetes$x, diabetes$y);
	    ITO <- MyLarsGiveMeClose(HitLars2Ob, IntClose =  5, MinMin = .00001);    


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
