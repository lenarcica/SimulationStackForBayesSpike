\name{CoordinateDescent}
\alias{CoordinateDescent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Coordinate Descent Lasso Operation}
\description{
   Wrapper for C++ code to fit Coordinate Descent Lasso Regression
}
\usage{
CoordinateDescent(xx = -1, yy = -1, XTX = -1, XTY = -1, NLen = -1,
   TotalLoops = -999, MaxEpsilon = -999, OnBeta = -1, OnGamma = -1, OnLambda = -1,
   RecordBetasFlag = FALSE, OnGammas = -999, InitKKs = -5,
   NumCDOConv = -999, CDOEpsilon = -999)
}
\arguments{
  \item{xx}{ covariate matrix, length NLen times klen }
  \item{yy}{ response vector, length NLen }
  \item{XTX}{ covariate information matrix transpose of xx times xx, length \eqn{\kappa \times \kappa}{} or klen*klen }
  \item{XTY}{ correlations vector, transpose of xx times yy, length klen }
  \item{NLen}{ length of vector yy if not supplied }
  \item{TotalLoops}{ Max number of times to run coordinate descent 
                                     loop on all covariates }
  \item{MaxEpsilon}{ Total criterion per loop to define convergence }
  \item{OnBeta}{ Initial vector to fit }
  \item{OnGamma}{ Single gamma penalty for all covariates }
  \item{OnLambda}{ Like gamma, only twice as large }
  \item{RecordBetasFlag}{ Flag that defines whether Coordinate Descent 
                               saves more response info}
  \item{OnGammas}{Different gamma parameters for all j in 1...klen covariates }
  \item{InitKKs}{ Toggles whether to do dynamic memory coordinate descent, 
                        faster for few active covariates
                  InitKKs should be a positive number.  Coordinate Descent will
                  Load initally only space for InitKKs vectors of XTX matrix 
                  into memory.   As vectors enter the model. The vector
                  t(X) \%*\% X_,j must be calculated.}
  \item{NumCDOConv}{Same as TotalLoops }
  \item{CDOEpsilon}{Same as MaxEpsilon }
}
\details{
  Coordinate Descent Lasso algorithm, uses C++ code
  
  Coordinate Descent is an iterative algorithm that loops over coordinates j in 
     1 ... kLen.  

  Friedman recommended a TotalLoops of order 30 to generate convergence in the
     kLen = 10, NLen = 350 point Lars Diabetes dataset.
     
  "OnGamma" is different from the Lambda parameter because the original paper
     defines a problem
   "Beta_Lasso" = ArgMin{  1/2 * sum_i( (Y_i - X_i * Beta)^2 ) +
                             gamma * sum_{j} abs(Beta_j) 
                        }
  This Gamma is arguably more intuitive than the LARS Lambda, and can be
    interpretted as Posterior behavior when the noise level is 1 and the 
    Laplacian coefficient is Gamma.  
}
\value{
  Returns an object with many quantities  
    \item{ReturnBetas}{Result vector}
    \item{FinalBeta}{Also final result vector}
    \item{XTX}{t(xx) \%*\% xx matrix}
    \item{XTY}{t(xx) \%*\% yy vector}
    \item{OnRecordBetas}{Matrix delineating output from all steps of ccordinate
           descent, useful for checking that TotalRuns/NumCDOConv parameter
           is large enough}
}
\references{ Friedman et al (2006)  "Pathwise Coordinate estimation" }
\author{ Alan Lenarcic }
\note{ "OnGammas" can be a vector with weights for each coefficient
}
\seealso{ CoordinateDescentPrototype }
\examples{
  ###### The Lars Analysis Format for the Paper
library(lars)
data(diabetes)
StlambdaA = 10; StlambdaD = StlambdaA;
XX <- StandardizeXX(diabetes$x);
YY <- StandardizeXX(diabetes$y);
MyLars <- lars(XX,YY);
## Member E is out
T1 <- proc.time();
CDO <- CoordinateDescent(xx= XX, yy = YY, OnGamma = 5, 
        NumCDOConv = 200, CDOEpsilon = .00001) 
CDO$FinalBeta;

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
