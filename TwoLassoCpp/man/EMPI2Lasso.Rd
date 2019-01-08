\name{EMPI2Lasso}
\alias{EMPI2Lasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ EMPI2Lasso Algorithm }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EMPI2Lasso(xxs = -1, yys = -1, NLen = -1, StlambdaD = -999, StlambdaA = -999, lambdaDMultC = 2^(2/8), lambdaAMultC = 2^(-2/8), lambdaAmultStop = 20, TotalRuns = 20, NumEMConv = 4, MultEMCons = 0.99, BetOlds = -999, StandardFlag = 0, RatWant = 0.2, pPpiuseSt = -999, sigmaNoiseSqSt = -999, m1 = 1, m2 = 1, NumEConv = 6, SigmaEta = -1, SigmaBarSq = -999, DConfidenceInt = -1, NumCDOConv = 50, CDOEpsilon = 1e-06, LambdaAK = -999, LambdaDK = -999, OrderSeq = -999, XTX = -1, XTY = -1, InverseGammaConstant = 1, WLSWeights = -1, TDFNu = -1, PrintFlag = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xxs}{ Covariate matrix, dimensions NLen * klen }
  \item{yys}{ Response vector, length NLen * 1 }
  \item{NLen}{ Length of yys vector if not supplied }
  \item{StlambdaD}{ Initial \eqn{\lambda_{\mbox{\tiny{$\mathcal{A}$}}}}{Lambda_A} value for input }
  \item{StlambdaA}{ Initial \eqn{\lambda_{\mbox{\tiny{$\mathcal{D}$}}}}{Lambda D} value for input }
 \item{lambdaDMultC}{ Factor to multiply \eqn{\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^{[\tau]}}{Lambda_D^[tau]} 
                 by to get \eqn{\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^{[\tau+1]}}{LambdaD^[tau+1] } }
  \item{lambdaAMultC}{ Factor to multiply 
     \eqn{ \lambda_{\mbox{\tiny{$\mathcal{A}$}}}^{[\tau]} }{Lambda_A^[tau]}
          by to get 
     \eqn{\lambda_{\mbox{\tiny{$\mathcal{A}$}}}^{[\tau+1]}}{Lambda_A^[tau+1]} }
  \item{lambdaAmultStop}{ Total number of iterations to decrease 
        \eqn{ \lambda_{\mbox{\tiny{$\mathcal{A}$}}}^{[\tau]} }{Lambda_A^[tau]} }
  \item{TotalRuns}{ Total number of iterations to increase 
      \eqn{\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^{[\tau]}}{Lambda_D^[tau]}}
 \item{NumEMConv}{ Number of EM converge loops, usually 5-8 suffices }
  \item{MultEMCons}{ Allows one to reduce amount of EM loops if requested }
 \item{BetOlds}{ Input Beta vector if desired }
  \item{StandardFlag}{ Standardize Covariates and response if desired }
  \item{RatWant}{ Not used anymore }
  \item{pPpiuseSt}{ Start value for seeking \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}}}{pi_A} }
  \item{sigmaNoiseSqSt}{ Start Value for seeking \eqn{\sigma^2}{sigma^2} }
  \item{m1}{ Together m1 and m2 are Beta prior constants for \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}}}{pi_A}  }
  \item{m2}{ Together m1 and m2 are Beta prior constants for \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}}}{pi_A} }
  \item{NumEConv}{ Old number for NumEMConv }
  \item{SigmaEta}{ Together SigmaEta and SigmaBarSq are Gamma distribution prior constants for \eqn{\sigma^2}{sigma^2} }
  \item{SigmaBarSq}{ Together SigmaEta and SigmaBarSq are Gamma distribution prior constants for \eqn{\sigma^2}{sigma^2} }
  \item{DConfidenceInt}{ Used in Confidence Interval Estimation but not operable. }
  \item{NumCDOConv}{ Maximum number of times to run Coordinate Descent }
  \item{CDOEpsilon}{ Sum of deviations to quit Coordinate Descent Algorithms }
 \item{LambdaDK}{ Vector of \eqn{\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^[\tau]}{Lambda_D^[tau]} values if that is more practical }
  \item{LambdaAK}{  Vector of \eqn{\lambda_{\mbox{\tiny{$\mathcal{A}$}}}^[\tau]}{Lambda_A^[tau]} values if that is more practical }
  \item{OrderSeq}{ Sequence of integers length TotalRuns (+4,-4,5,-5,5...) For each value of 
     \eqn{\left(\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^[\tau], 
        \lambda_{\mbox{\tiny{$\mathcal{A}$}}}^[\tau] \right)}{(Lambda_D^[t],Lambda_A^[t])} this decides whether to
         do E step first or M step (sign) and also how many EM steps per value (magnitude) }
 \item{XTX}{ transpose of xx times xx if if desired }
  \item{XTY}{ transpose of xxs times yys if desired to pass less data }
\item{InverseGammaConstant}{ Changes memory management within function when it inverts the gamma weights }
  \item{WLSWeights}{ If a positive vector of length yys, weights the data for Weighted Least Squares }
  \item{TDFNu}{ If a positive value, then ammends EM to weight Y data accoridng to t-distribution Noise of given coefficient.  Warning
   that in general, low degree of freedom T distributions produce multi-modal solutions, even in non-selection problems.  Limit-Lasso
   approach might be most helpful in generating large mode.}
}
\details{
   Estimates EM2Lasso fit while giving uncertain beta prior distribution for 
     "pi_A" parameter
}
\value{
    Returns TwoLasso object  
     \item{ReturnBetasAll}{Gives the final betas for all input values of LambdaA/LambdaD pairs}
     \item{BBHatsAll}{gives theestimates for expected B indicator values.}
     \item{ReturnBetas}{Final Result output from Limit-Lasso final LambdaA/K values}
     \item{LambdAK}{Used Sequence of 
        \eqn{\lambda_{\mbox{\tiny{$\mathcal{A}$}}}^{[\tau]}}{Lambda_A^[tau]} 
        parameters}
     \item{LambdaDK}{ Used Sequence of
           \eqn{\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^{[\tau]}}{Lambda_D^[tau]} 
        parameters
     }
     \item{PiRecVec}{Used Sequence of \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}}}{pi_A}
             parameters for Limit Lasso Loops}
     \item{USReturnBetasAll}{If StandardFlag = 1, then this is the 
           "Before unstandardization" vales of the beta coefficients}
     \item{sdXXs}{Standard deviations of XXs before standardization}
     \item{sdYYs}{Standard deviations of YYs before standardization}
     \item{USsigmaNoiseSq}{sigmaNoise used, assuming before unstandardization}
     \item{FailureFlag}{If something has gone wrong, FailureFlag will be negative}
}
\references{ Alan Lenarcic's 2009 Harvard Thesis Chapter 2

   "TwoLasso for Bayes Adapted Selection", Lenarcic, Journal of
   Computational and Graphical Statistics, Submitted 2009. }
\author{ Alan Lenarcic }
\note{   See EM2Lasso for ideas
}
\seealso{ EM2Lasso, EMRIDGE }
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
EMPI2LassoToy1 <- EMPI2Lasso(xxs= XX, yys = YY, pPpiuseSt = .5, 
         sigmaNoiseSqSt = .5,  SigmaEta = 40, SigmaBarSq = .5,
           RatWant = .2, 
          StlambdaD = StlambdaD, StlambdaA = StlambdaA, lambdaDMultC = 3^(1/8), 
          lambdaAMultC = 1/3^(1/8), lambdaAmultStop = 12, TotalRuns = 12, 
          NumEMConv = 7, MultEMCons = .99, BetOlds = -999, StandardFlag = 1,
          m1 = 10, m2 = 10) 
  t(EMPI2LassoToy1$ReturnBetasAll[, ]);
  t(EMPI2LassoToy1$BBHatsAll[, ]);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
