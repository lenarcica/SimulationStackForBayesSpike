\name{SimMeStuff}
\alias{SimMeStuff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ SimMeStuff }
\description{
  Simulation Code for Estimator Comparison
}
\usage{
SimMeStuff(NN = -999, NR = -999, Noise = -999, BRealVecs = -999, ppCov = -999, sigmaFactCov = -999)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NN}{ Number of Sample Points}
  \item{NR}{ Number of Covariates }
  \item{Noise}{ Noise Level "Sigma", not SigmaSquared }
  \item{BRealVecs}{ Real values for true nonzero Beta Coefficients }
  \item{ppCov}{ Probability of covariance between two covariates }
  \item{sigmaFactCov}{ Distribution of covariance between two covariates }
}
\details{
   Simulates Some Data sets for Estimator Comparison
}
\value{
  A lot of stuff
}
\references{ "TwoLasso for Bayes Adapted Selection", Lenarcic, Journal of
   Computational and Graphical Statistics, Submitted 2009. }
\author{ Alan Lenarcic }
\note{  Nothing Interesting
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
SimMeStuff(NN=10, NR = 6, Noise = .5, BRealVecs=c(4,3,-2.5,1),
    ppCov = .25, sigmaFactCov = .3);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
