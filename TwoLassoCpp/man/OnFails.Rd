\name{OnFails}
\alias{OnFails}
\docType{data}
\title{ Number counting simulation fit failures }
\description{
  An integer, 0 if there have been no SimForMe2 failures.
}
\usage{OnFails}
\format{
  The format is:
 num 0
}
\details{
   Number of Total Fails so far
}
\source{ Alan Lenarcic
}
\seealso{ FailToFitSim, FailToFitEMOb }
\references{ "TwoLasso for Bayes Adapted Selection", Lenarcic, Journal of
   Computational and Graphical Statistics, Submitted 2009.
}
\examples{
#print(OnFails);
## maybe str(OnFails) ; plot(OnFails) ...
}
\keyword{datasets}
