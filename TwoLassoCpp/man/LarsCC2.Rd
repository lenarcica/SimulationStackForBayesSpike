\name{LarsCC2}
\alias{LarsCC2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ LARS algorithm coded in CC }
\description{
  LARS regression Algorithm for Lasso
  
  Performs Weighted LARS as well as unwiehged lars.
  
  If StandardFlag is 0 then there will be no standardization before algorithm.
}
\usage{
LarsCC2(xxs, yys, BetStart = -1, lambda = 0, BetOlds = -1, StandardFlag = 0, Weights = NULL, GFlag = -1, PrintGFlag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xxs}{ x covariate matrix }
  \item{yys}{ y response vector}
  \item{BetStart}{ Initiate Beta vector }
  \item{lambda}{ Final smallest \eqn{\lambda}{lambda} to stop LARS algorithm on }
  \item{BetOlds}{ Same as BetStart }
  \item{StandardFlag}{ Set to 1 to standardize X and Y before algorithm }
  \item{Weights}{ Weight vetor to give to covariates }
  \item{GFlag}{ Decides whether to create "covariance" matrix G = t(X)*X, 
     which toggles whether
      LARS is an order n *k_a^3 algorithm or n*k*k + k_a^3 }
  \item{PrintGFlag}{ Set flag to a number to see printing  during algorithm}
}
\details{
   Runs LARS algorithm as coded in C++.  Stops when reaching lambda for stopping point.
}
\value{
  \item{ReturnBetas}{Final Values of Beta coefficients for fixed inputted lambda parameter}
  \item{InitBetas}{First Betas given to start algorithm (should be all zeroes if you want 
     right answer)}
  \item{BetasKeep}{matrix of Beta coefficients along the LARS path starting from 0}
  \item{KeepJs}{Count of active coefficients along each step of LARS path}
  \item{Penalty}{Matrix keeping measure of current sum of squares and current
           sum of absolute values that mark the optimal Lasso "posterior"
            " sum( (Y_i-X_i * beta)^2 ) + sum(lambda * abs(beta_j)) "
        }
  \item{OnBetas}{Final Beta vector settled upon}
  \item{lambda}{\eqn{\lambda}{lambda} parameter that algorithm stopped upon}
  \item{USReturnBetas}{When StandardFlag = 1 then X and Y were standardized.
      Before the ReturnBetas output from LARS is unstandardized (i.e. multiplied by
       sdYYs/sdXXs), this is LARS beta coefficients}
  \item{USWReturnBetas}{When Weights are nonzero then. We convert Betas to an
     even weighted scale.  This is ReturnBetas 
     output from LARS algorithm before reweighting and 
     unstandardization.
  }
  \item{USReturnBetasAll}{When StandardFlag = 1 then X and Y were standardized.
      Before the ReturnBetasAll output
      (All steps up to \eqn{\lambda}{lambda}) from LARS is unstandardized (i.e. multiplied by
       sdYYs/sdXXs), this is LARS beta coefficients}
  \item{USWReturnBetasAll}{When Weights are nonzero then. We convert Betas to an
     even weighted scale.  This is ReturnBetasAll 
     (All steps up to \eqn{lambda}{lambda}) 
     output from LARS algorithm before reweighting and 
     unstandardization.
  \item{Lambdas}{Should be the Lambdas that cause branching points up to the final
     \eqn{\lambda}{lambda}.  However, those seem to be difficult to calculate, and 
     I do not calculate for true \eqn{\lambda}{lambda} when there are no non-zero coefficients.}
  }
}
\references{ Efron LARS algorithm 2004 }
\author{Alan Lenarcic }
\note{
}
\seealso{ LarsPrototype, CoordinateDescent }
\examples{
#### Lars Scratch
###
### Our Goal eventually is algorithm that minimizes
### sum( Y - XB )^2 + v_1 ||B||_1 + v_2|| B- B_0||_1
NN <- 20;
Noise <- .1
BReal <- matrix( c(1, 0, 1.5,-.5, 0), 5,1);
LCovX <- rbind( c(1,0,0,0,0),
                c(.2, sqrt(1-.2^2),0,0,0),
                c(.1,.3, sqrt( 1- .1^2-.3^2), 0,0),
                c(-.25,.1,.1, sqrt( 1- .25^2 - .1^2 -.1^2 ),0),
                c(-.3, .05, .05, .2, sqrt( 1- .3^2 -.05^2-.05^2 - .2^2) )
              )
CovX <- LCovX \%*\% t(LCovX)
XX <- LCovX  \%*\% matrix( rnorm(NN * length(CovX[1,]),0,1), length(CovX[1,]), NN);
XX <- t(XX);
YY <-    rnorm(NN,0,Noise) + XX \%*\% BReal;


StandardizeXX <- function(XX) {
  if ( length(dim(XX)) == 2) {
        XX <- t( ( t(XX) - meanC(XX)) / sd(XX) );
  }  else {
        XX <-  (XX - mean(XX)) / sd(XX);
  }
  return(XX);
}
yys <- StandardizeXX(YY);
OnBetas <- BReal * c(1,1,0,0,0);
xxs <- StandardizeXX(XX);
OnActives <-  c(1,1,0,0,0);
BetOlds <- c(0,0,0,0,0);

 yys <- YY - mean(YY);
 yys = matrix(yys,length(yys), 1);
 xxs <- XX;
 xxs = matrix(xxs, length(xxs[,1]), length(xxs[1,]) );
lambda <- .01;
BetStart = 0 * OnBetas; lambda = 1; BetOlds = -1;
StandardFlag = 0; Weights = NULL; GFlag = -1; PrintGFlag = FALSE;
     Blasso <- LarsCC2(xxs=xxs,yys=yys,  BetStart=OnBetas* 0, lambda = lambda, BetOlds = -1);
                  Beta <- Blasso$OnBetas;
    ## Alars <- lars(x=xxs,y=as.vector(yys))
          Blasso <- LarsCC2(xxs=xxs,yys=yys,  OnBetas* 0, lambda = lambda, BetOlds = -1, StandardFlag = 1);
                  Beta <- Blasso$OnBetas;
    ## Alars <- lars(XX,as.vector(YY))
    ##  plot(Alars)
  Blasso <- LarsCC2(xxs,yys,  OnBetas* 0, lambda = lambda,  BetOlds = -1);
                  Beta <- Blasso$OnBetas;     
   ##Alasso2 <- l1ce(yys~xxs-1, bound = sum(Blasso$OnBetas) );
   Classo <- LarsPrototype(xxs=xxs,yys=yys, BetStart = OnBetas * 0, 
         lambda = lambda, BetOlds = -1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
