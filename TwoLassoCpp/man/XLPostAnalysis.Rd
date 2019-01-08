\name{XLPostAnalysis}
\alias{XLPostAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Marginal Median Posterior Analysis }
\description{
   Given a Sparse REML estimate, such as Limit-Lasso input, this calculates
    the "Marginal Median Estimate"
}
\usage{
XLPostAnalysis(BetaProposed, xxs = -1, yys = -1, NLen, ppiuse, sigmaNoiseSq, tausqA, XTX = -1, XTY = -1, YtY = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BetaProposed}{ ~~Describe \code{BetaProposed} here~~ }
  \item{xxs}{ ~~Describe \code{xxs} here~~ }
  \item{yys}{ ~~Describe \code{yys} here~~ }
  \item{NLen}{ ~~Describe \code{NLen} here~~ }
  \item{ppiuse}{ ~~Describe \code{ppiuse} here~~ }
  \item{sigmaNoiseSq}{ ~~Describe \code{sigmaNoiseSq} here~~ }
  \item{tausqA}{ ~~Describe \code{tausqA} here~~ }
  \item{XTX}{ ~~Describe \code{XTX} here~~ }
  \item{XTY}{ ~~Describe \code{XTY} here~~ }
  \item{YtY}{ ~~Describe \code{YtY} here~~ }
}
\details{
  Marginal density uses the point estimate for other "klen-1" coordinates.
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ TwoLasso for Bayes Adapted Selection, Lenarcic, JCGS 2009 }
\author{ Alan Lenarcic }
\note{ 
}
\seealso{ SimForMe2 }
\examples{
library(lars)
data(diabetes)
StlambdaA = 10; StlambdaD = StlambdaA;
XX <- StandardizeXX(diabetes$x);
YY <- StandardizeXX(diabetes$y);
MyLars <- lars(XX,YY);
## Member E is out
T1 <- proc.time();
EM2LassoToy1 <- EM2Lasso(xxs= XX, yys = YY, ppiuse = .5, 
         sigmaNoiseSq = .5, 
           RatWant = .2, 
          StlambdaD = StlambdaD, StlambdaA = StlambdaA, lambdaDMultC = 3^(1/8), 
          lambdaAMultC = 1/3^(1/8), lambdaAmultStop = 12, TotalRuns = 12, 
          NumEMConv = 7, MultEMCons = .99, BetOlds = -999, StandardFlag = 1) 
  t(EM2LassoToy1$ReturnBetasAll[, ]);
  t(EM2LassoToy1$BBHatsAll[, ]);
  XLPostAnalysis(EM2LassoToy1$ReturnBetas, xxs = XX, yys = YY, 
        NLen = length(YY), ppiuse = .5, sigmaNoiseSq = .5,
        tausqA = 5);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
