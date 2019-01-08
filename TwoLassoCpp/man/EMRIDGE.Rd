\name{EMRIDGE}
\alias{EMRIDGE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EMRIDGE(xxs, yys, TotalRuns = 128, NumEMConv = 4, ppiuse = 0.5, tauDsqStart = -999, tauAsqStart = -999, tauDMult = -999, tauAMult = -999, tauAMultStop = 0, BStartBetas = -1, sigmaNoiseSq = 0.5, SigmaSqNoise = -9, EFlag = FALSE, logp = -1, log1mp = -1, SigmaNoiseInv = -1, NOPRINT = TRUE, n1 = -9, n2 = -9, NumCDOConv = 40, CDOEpsilon = 1e-05, TauSqASeq = -1, TauSqDSeq = -1, OrderSeq = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xxs}{ Covariate matrix, dimensions NLen * klen }
  \item{yys}{ Response vector, length NLen * 1 }
  \item{TotalRuns}{ Total number of iterations to increase \eqn{\tau_{\mbox{\tiny{$\mathcal{D}$}}}^2[t]}{tau_D^2[t]}}
  \item{NumEMConv}{ Number of EM converge loops, usually 5-8 suffices }
  \item{ppiuse}{ \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}}}{pi_A} probability  }
  \item{tauDsqStart}{ Initial \eqn{\tau_{\mbox{\tiny{$\mathcal{D}$}}}^2}{tau_D^2} value for input }
  \item{tauAsqStart}{ Initial \eqn{\tau_{\mbox{\tiny{$\mathcal{A}$}}}^2}{tau_A^2} value for input }
  \item{tauDMult}{ Multiplier to turn \eqn{{\tau_{\mbox{\tiny{$\mathcal{D}$}}}^2}^{[t]}}{tau_D^2^[t]} 
      into \eqn{{\tau_{\mbox{\tiny{$\mathcal{D}$}}}^2}^{[t]}}{tau_D^2^[t+1]} }
  \item{tauAMult}{ Multiplier to turn \eqn{{\tau_{\mbox{\tiny{$\mathcal{A}$}}}^2}^{[t]}}{tau_A^2^[t]} 
      into \eqn{{\tau_{\mbox{\tiny{$\mathcal{A}$}}}^2}^{[t]}}{tau_A^2^[t+1]} }
  \item{tauAMultStop}{ Last iteration to up \eqn{{\tau_{\mbox{\tiny{$\mathcal{A}$}}}^2}^{[t]}}{tau_A^2^[t]}  }
  \item{BStartBetas}{ Starting vector for \eqn{\boldsymbol{\beta}}{Beta} Vector }
  \item{sigmaNoiseSq}{ Value for \eqn{\sigma^2}{sigma^2} }
  \item{SigmaSqNoise}{ Value for \eqn{\sigma^2}{sigma^2}  }
  \item{EFlag}{ A flag for returning E step information  }
  \item{logp}{ An input of \eqn{\log_{e} \left\{ \pi_{\mbox{\tiny{$\mathcal{A}$}}} \right\}}{log pi_A}, not needed }
  \item{log1mp}{ An input of \eqn{\log_{e} \left\{ 1-\pi_{\mbox{\tiny{$\mathcal{A}$}}} \right\}}{log (1-pi_A)}, not needed }
  \item{SigmaNoiseInv}{ An input of \eqn{1/\sigma^2}{1/sigma^2}, not needed. }
  \item{NOPRINT}{  Flag to not print information. }
  \item{n1}{ Bayesian weight, n1 and n2 are Beta prior weights instead of \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}}}{pi_A} }
  \item{n2}{ Bayesian weight, n1 and n2 are Beta prior weights instead of \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}}}{pi_A} }
 \item{NumCDOConv}{ Maximum number of times to run Coordinate Descent, if desired instead of matrix inversion }
  \item{CDOEpsilon}{ Sum of deviations to quit Coordinate Descent Algorithms }
  \item{TauSqASeq}{ A sequence of input \eqn{\tau_{\mbox{\tiny{$\mathcal{A}$}}}^2[t]}{tau_A^2[t]} values if desired }
  \item{TauSqDSeq}{ A sequence of input \eqn{\tau_{\mbox{\tiny{$\mathcal{D}$}}}^2[t]}{tau_D^2[t]} values if desired }
  \item{OrderSeq}{ Sequence of integers length TotalRuns (+4,-4,5,-5,5...) For each value of 
     \eqn{\left({\tau_{\mbox{\tiny{$\mathcal{D}$}}}^2}^[t], 
        {\tau_{\mbox{\tiny{$\mathcal{A}$}}}^2}^[t] \right)}{(tau_D^2^[t],tau_A^2^[t])} this decides whether to
         do E step first or M step (sign) and also how many EM steps per value (magnitude) }
}
\details{
   Treat like a Limit-Lasso regression, selection only occurs if
     TauSqASeq and TauSqDSeq are shrunk far.
}
\value{
   Returns EMLARS object, ReturnBetas is Return vector.  BBHatsFinal is expected
     B values.
}
\references{ Alan Lenarcic's Harvard Thesis 2009 Chapter 2

   "TwoLasso for Bayes Adapted Selection", Lenarcic, Journal of
   Computational and Graphical Statistics, Submitted 2009.}
\author{ Alan Lenarcic }
\note{  Not Really
}
\seealso{ EMLARS, CoordinateDescent, EMPILARS }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

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

EM2RIDGEToy1 <- EMRIDGE(xxs= XX, yys = YY, ppiuse = .5, 
         sigmaNoiseSq = .5, 
         TauSqASeq = .5 * 2^((0:64)/8),
         TauSqDSeq = .25 * 2^(-(0:64)/8),
         OrderSeq = rep(6, 65),
         NumEMConv = 7,
         CDOEpsilon = 1e-05, NumCDOConv = 40, NOPRINT = FALSE) 
  (EM2RIDGEToy1$ReturnBetasAllEM[,1,c(1,(1:8)*8) ]);
  (EM2RIDGEToy1$BBHatsAllEM[,1,c(1,(1:8)*8) ]);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
