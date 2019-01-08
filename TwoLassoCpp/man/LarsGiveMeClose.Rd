\name{LarsGiveMeClose}
\alias{LarsGiveMeClose}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LarsGiveMeClose(LarsOb, IntClose, MinMin = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LarsOb}{ Object outputted by lars program.  }
  \item{IntClose}{ Number of coefficients to seek for active }
  \item{MinMin}{ Minimum value for Lambda parameter }
}
\details{
  Given a Lars Object (outed from program lars preferentially), 
    seeks the first branching step such that "IntClose" coefficients are active.
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
\note{
}
\seealso{ LarsCC2, lars, CoordinateDescent }
\examples{
      library(lars); data(diabetes);
	    HitLars2Ob <- lars(diabetes$x, diabetes$y);
	    ITO <- LarsGiveMeClose(HitLars2Ob, IntClose =  5, MinMin = .00001);    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
