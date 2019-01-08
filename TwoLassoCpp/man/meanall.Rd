\name{meanall}
\alias{meanall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mean of All Columns }
\description{
  Also just colMeans
}
\usage{
meanall(MMT)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MMT}{ A Matrix of numbers }
}
\details{
  Same as colMeans
}
\value{
Returns vector with means of each column of MMT
}
\references{ "TwoLasso for Bayes Adapted Selection", Lenarcic,
  Journal of Computational and Graphical Statistics, Submitted 2009. }
\author{ Alan Lenarcic }
\note{  Same as colMeans
}
\seealso{ colMeans, meanC }
\examples{
meanall(matrix(1:25,5,5));
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ arith }% __ONLY ONE__ keyword per line
