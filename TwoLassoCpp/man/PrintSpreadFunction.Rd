\name{PrintSpreadFunction}
\alias{PrintSpreadFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Coordinate Descent Lasso Operation}
\description{
  Prints Trace of BBhat estimates for EM2Lasso Regression
}
\usage{
PrintSpreadFunction(EMObject, ITFlag = TRUE, weightFlag = FALSE, 
    USFlag = FALSE, logFlag = FALSE)
}
\arguments{
  \item{EMObject}{ Output from EM2Lasso, EMRIDGE, or EM2PILasso}
}
\details{
  Other Flags effect the header to table
}
\value{
   Nothing
}
\references{ Lenarcic 2009 Thesis }
\author{ Alan Lenarcic }
\note{ Can tell whether one did standardized regression or no. 
}
\seealso{ EM2Lasso }
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
T2 <- proc.time();
### Member E stays
  EM2LassoToy1 <- EM2Lasso(xxs= XX, yys = YY, ppiuse = .6, 
         sigmaNoiseSq = .5, 
           RatWant = .2, 
          StlambdaD = StlambdaD, StlambdaA = StlambdaA, lambdaDMultC = 3^(1/8), 
          lambdaAMultC = 1/3^(1/8), lambdaAmultStop = 12, TotalRuns = 24, 
          NumEMConv = 7, MultEMCons = .99, BetOlds = -999, StandardFlag = 1) 
  t(EM2LassoToy1$ReturnBetasAll[, ]);
  t(EM2LassoToy1$BBHatsAll[, ]);
      PrintSpreadFunction(EM2LassoToy1)
  
#### Now Member E Stays
     StlambdaA = 10; StlambdaD = StlambdaA;
     XX <- StandardizeXX(diabetes$x);
     YY <- StandardizeXX(diabetes$y);
     MyLars <- lars(XX,YY);
     ## Member E is out
     T1 <- proc.time();
StlambdaD = 25; StlambdaA = 25; ppiuse = .6
     EM2LassoToy1 <- EM2Lasso(xxs= XX, yys = YY, ppiuse = ppiuse, 
              sigmaNoiseSq = .25, 
                RatWant = .2, 
               StlambdaD = StlambdaD, StlambdaA = StlambdaA, lambdaDMultC = 3^(.1/8), 
               lambdaAMultC = 1/3^(.1/8), lambdaAmultStop = 50, TotalRuns = 400, 
               NumEMConv = 7, MultEMCons = .99, BetOlds = -999, StandardFlag = 1, 
              OrderSeq = rep(-8, 400) ) 
       t(EM2LassoToy1$ReturnBetasAll[c(1,1:10 * 40), ]);
       t(EM2LassoToy1$BBHatsAll[c(1,1:10 * 40), ]); 
      PrintSpreadFunction(EM2LassoToy1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ hplot }% __ONLY ONE__ keyword per line
