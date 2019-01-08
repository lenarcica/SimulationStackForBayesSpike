\name{LarsPrototype}
\alias{LarsPrototype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ R Coded Lars Algorithm }
\description{
  R Protottype Function runs LARS Algorithm
}
\usage{
LarsPrototype(xxs, yys, BetStart = -1, lambda = 0, BetOlds = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xxs}{ x covariate matrix }
  \item{yys}{ y response vector}
  \item{BetStart}{ Initiate Beta vector -- Not Very Useful for LARS, though good for testing }
  \item{lambda}{ Final smallest lambda to stop LARS algorithm on }
  \item{BetOlds}{ Initial BetaVector, same as BetStart }
}
\details{
  Runs LARS Algorithm starting from BetStart (or zero default) until reaching Lambda
}
\value{
 Returns EMLARS object.  ReturnBetas is return vector
}
\references{  }
\author{ Alan Lenarcic }
\note{    See LarsCC2, hardcoded version
}
\seealso{ LarsCC2 }
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
CovX <- LCovX %*% t(LCovX)
XX <- LCovX  %*% matrix( rnorm(NN * length(CovX[1,]),0,1), length(CovX[1,]), NN);
XX <- t(XX);
YY <-    XX %*% BReal + rnorm(NN, 0, Noise );

yys <- StandardizeXX(YY);
OnBetas <- BReal * c(1,1,0,0,0);
xxs <- StandardizeXX(XX);
OnActives <-  c(1,1,0,0,0);
BetOlds <- c(0,0,0,0,0);

 yys <- YY - mean(YY);
 xxs <- XX;
lambda <- .01;
     Blasso <- LarsCC2(xxs,yys,  OnBetas* 0, lambda = lambda, BetOlds = -1);
                  Beta <- Blasso$OnBetas;
     ##Alars <- lars(xxs,yys)
          Blasso <- LarsCC2(XX,YY,  OnBetas* 0, lambda = lambda,  BetOlds = -1, StandardFlag = 1);
                  Beta <- Blasso$OnBetas;
     ##Alars <- lars(XX,YY)
     ##plot(Alars)
  Blasso <- LarsCC2(xxs,yys,  OnBetas* 0, lambda = lambda, BetOlds = -1);
                  Beta <- Blasso$OnBetas;     
   ##Alasso2 <- l1ce(yys~xxs-1, bound = sum(Blasso$OnBetas) );
   Classo <- LarsPrototype(xxs,yys, BetStart = OnBetas * 0, lambda = lambda, 
       BetOlds = -1)
   }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
