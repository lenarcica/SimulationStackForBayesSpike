\name{TwoLassoCpp-package}
\alias{TwoLassoCpp-package}
\alias{TwoLassoCpp-package}
\alias{TwoLasso}
\docType{package}
\title{
TwoLassoCpp Package

A \code{\link[Rcpp]{Rcpp}} package using Two Lasso objects.  Uses Coordinate descent algorithms
to dynamically solve for 2 prior selection comparison.  Essentially a fixed point approximation
for methods used in package \code{\link[BayesSpike]{BayesSpike}} which is a MCMC integrative package.
 
 See \code{\link{TwoLassoCpp}()} for info on the main Two Lasso function to call.
 
 Code exists in this package on advice to other possible approaches and other attempts at implementing TwoLasso
 with Coordinate Descent or other approaches.
}
\description{
Implementation of TwoLasso, LimitLasso and related algorithms
}
\details{
\tabular{ll}{
Package: \tab TwoLasso\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-05-07\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
Important Functions:
  EM2Lasso, EMPI2Lasso
}
\author{
Alan Lenarcic

Maintainer: Who to complain to <alenarc@stat.harvard.edu>

}
\references{
  TwoLasso for Bayes Adapted Selection
}
\keyword{ methods }
\keyword{ multivariate }
\seealso{ \code{\link{TwoLassoCpp}}, 
  \code{\link{TwoLassoSEXP}}, \code{\link{TwoLassoRegression}},
  \code{\link{TwoLassoCrossValidate}}, \code{\link{EM2Lasso}}, 
  \code{\link{EMPI2Lasso}}, \code{\link{EMRIDGE}}, 
  \code{\link{XLPostOperations}}
}
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
}