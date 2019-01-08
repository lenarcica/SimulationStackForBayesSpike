\name{GLMG2Lasso}
\alias{GLMG2Lasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ GLMG2Lasso Algorithm }
\description{
  On Logistic Data,  Runs Two-Lasso algorithm and Limit-Lasso using EM and LARS/Coordinate Descent to do Lasso optimization
}
\usage{
GLMG2Lasso(X, y01, StartBeta = -999, StartBeta0 = 0, piA = 0.5, 
    sigmaSq = 2, LambdaAK = -999, LambdaDK = -999, OrderSeq = -999, 
    StLambdaA = 1, StLambdaD = 1, MultLambdaA = 0.98, MultLambdaD = 0.98, 
    TotalRuns = 100, MaxLambdaASteps = 50, RecordBetaFlag = TRUE, 
    InitKKs = 5, PrintFlag = 0, FixKa = -1) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ Covariate matrix, dimensions NLen * klen }
  \item{y01}{ 0 or 1 Response vector, length NLen * 1 }
  \item{StartBeta}{ Starting values of Beta to start with }
  \item{piA}{ \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}}}{pi_A} probability  }
  \item{sigmaSq}{ \eqn{\sigma^2}{sigma^2} input value, This sigma value acts as a smoother
     for parameter acceptance between logistic regression, default 2 kills that parameter.}
  \item{StLambdaD}{ Initial \eqn{\lambda_{\mbox{\tiny{$\mathcal{A}$}}}}{Lambda_A} value for input }
  \item{StLambdaA}{ Initial \eqn{\lambda_{\mbox{\tiny{$\mathcal{D}$}}}}{Lambda D} value for input }
  \item{MultLambdaA}{ Factor to multiply \eqn{\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^{[\tau]}}{Lambda_D^[tau]} 
                 by to get \eqn{\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^{[\tau+1]}}{LambdaD^[tau+1] } }
  \item{MultLambdaD}{ Factor to multiply 
     \eqn{\lambda_{\mbox{\tiny{$\mathcal{A}$}}}^{[\tau]}}{Lambda_A^[tau]} by to get 
     \eqn{\lambda_{\mbox{\tiny{$\mathcal{A}$}}}^{[\tau+1]}}{Lambda_A^[tau+1]} }
  \item{MaxLambdaASteps}{ Total number of iterations to decrease \eqn{\lambda_{\mbox{\tiny{$\mathcal{A}$}}}^{[\tau]}}{Lambda_A^[tau]} }
  \item{TotalRuns}{ Total number of iterations to increase \eqn{\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^{[\tau]}}{Lambda_D^[tau]}}
  \item{MaxCDOLoops}{ Number of EM converge loops, usually 5-8 suffices }
  \item{MultEMCons}{ Allows one to reduce amount of EM loops if requested }
  \item{StandardFlag}{ Standardize Covariates and response if desired }
  \item{LambdaDK}{ Vector of \eqn{\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^[\tau]}{Lambda_D^[tau]} values if that is more practical }
  \item{LambdaAK}{  Vector of \eqn{\lambda_{\mbox{\tiny{$\mathcal{A}$}}}^[\tau]}{Lambda_A^[tau]} values if that is more practical }
  \item{OrderSeq}{ Sequence of integers length TotalRuns (+4,-4,5,-5,5...) For each value of 
     \eqn{\left(\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^[\tau], 
        \lambda_{\mbox{\tiny{$\mathcal{A}$}}}^[\tau] \right)}{(Lambda_D^[t],Lambda_A^[t])} this decides whether to
         do E step first or M step (sign) and also how many EM steps per value (magnitude) }
  \item{NumCDOConv}{ Maximum number of times to run Coordinate Descent }
  \item{CDOEpsilon}{ Sum of deviations to quit Coordinate Descent Algorithms }
  \item{InverseGammaConstant}{ Changes memory management within function when it inverts the gamma weights }
  \item{FixKa}{ If a positive nonzero value, implements Fermi-Dirac Approximation Limit Lasso, finding best step with FixKa active coefficients }
  \item{WLSWeights}{ If a positive vector of length yys, weights the data for Weighted Least Squares }
  \item{TDFNu}{ If a positive value, then ammends EM to weight Y data accoridng to t-distribution Noise of given coefficient.  Warning
   that in general, low degree of freedom T distributions produce multi-modal solutions, even in non-selection problems.  Limit-Lasso
   approach might be most helpful in generating large mode.}
}
\details{
   The EM2Lasso algorithm for your Two-Lasso and Limit-Lasso needs.
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
             parameters for Limit Lasso Loops.  Note that if "FixKa" was a positive
             integer, then \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}}}{pi_A} will be forced
             to vary throughout the algorithm. }
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
\note{ 
    Limit-Lasso is obtained by giving LambdaAK/LambdaDK sequences,
     or alternatively, giving starting StLambdaA/StLambdaD and multipliers.
}
\seealso{ EMPILARS, EMRIDGE, CoordinateDescent, CoordinateDescentPrototype }
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
EM2LassoToy1 <- EM2Lasso(xxs= XX, yys = YY, ppiuse = .5, 
         sigmaNoiseSq = .5, 
           RatWant = .2, 
          StlambdaD = StlambdaD, StlambdaA = StlambdaA, lambdaDMultC = 3^(1/8), 
          lambdaAMultC = 1/3^(1/8), lambdaAmultStop = 12, TotalRuns = 12, 
          NumEMConv = 7, MultEMCons = .99, BetOlds = -999, StandardFlag = 1) 
  t(EM2LassoToy1$ReturnBetasAll[, ]);
  t(EM2LassoToy1$BBHatsAll[, ]);
     PrintSpreadFunction(EM2LassoToy1);
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
  
#### Now Member E Stays
     StlambdaA = 10; StlambdaD = StlambdaA;
     XX <- StandardizeXX(diabetes$x);
     YY <- StandardizeXX(diabetes$y);
     MyLars <- lars(XX,YY);
     ## Member E is out

StlambdaD = 25; StlambdaA = 25; ppiuse = .5
     EM2LassoToy1 <- EM2Lasso(xxs= XX, yys = YY, ppiuse = ppiuse, 
              sigmaNoiseSq = .25, 
                RatWant = .2, 
               StlambdaD = StlambdaD, StlambdaA = StlambdaA, lambdaDMultC = 3^(.1/8), 
               lambdaAMultC = 1/3^(.1/8), lambdaAmultStop = 50, TotalRuns = 400, 
               NumEMConv = 7, MultEMCons = .99, BetOlds = -999, StandardFlag = 1, 
              OrderSeq = rep(-8, 400) ) 
       t(EM2LassoToy1$ReturnBetasAll[c(1,1:10 * 40), ]);
       t(EM2LassoToy1$BBHatsAll[c(1,1:10 * 40), ]);
   PrintSpreadFunction(EM2LassoToy1);
#######################3
### PrintSpread Demonstration side by side
StlambdaD = .01; StlambdaA = .01; ppiuse = .95
     EM2LassoToy1a <- EM2Lasso(xxs= XX, yys = YY, ppiuse = ppiuse, 
              sigmaNoiseSq = .5, 
                RatWant = .2, 
               StlambdaD = StlambdaD, StlambdaA = StlambdaA, 
               lambdaDMultC = 3^(.25/8), 
               lambdaAMultC = 1/3^(.1/8), lambdaAmultStop = 50, TotalRuns = 250, 
               NumEMConv = 14, MultEMCons = .99, BetOlds = -999, StandardFlag = 1, 
              OrderSeq = rep(14, 450), CDOEpsilon = .00000000001, NumCDOConv = 4000) 
       t(EM2LassoToy1$ReturnBetasAll[c(1,1:10 * 80), ]);
       t(EM2LassoToy1$BBHatsAll[c(1,1:10 * 80), ]);
   PrintSpreadFunction(EM2LassoToy1a);
  summary(lm(YY~XX[,c(2,3,4,5,6,9)]));
StlambdaD = .5; StlambdaA = .5; ppiuse = .75
     EM2LassoToy1b <- EM2Lasso(xxs= XX, yys = YY, ppiuse = ppiuse, 
              sigmaNoiseSq = .5, 
                RatWant = .2, 
               StlambdaD = StlambdaD, StlambdaA = StlambdaA, lambdaDMultC = 3^(.1/8), 
               lambdaAMultC = 1/3^(.1/8), lambdaAmultStop = 50, TotalRuns = 250, 
               NumEMConv = 14, MultEMCons = .99, BetOlds = -999, StandardFlag = 1, 
              OrderSeq = rep(14, 450), CDOEpsilon = .00000000001, NumCDOConv = 4000) 
       t(EM2LassoToy1$ReturnBetasAll[c(1,1:10 * 80), ]);
       t(EM2LassoToy1$BBHatsAll[c(1,1:10 * 80), ]);
   PrintSpreadFunction(EM2LassoToy1b);
  summary(lm(YY~XX[,c(2,3,4,5,6,9)])); 
  matlayout <- rbind(
     c(1,1,1),
     c(3,2,5),
     c(4,2,6)
     )
nf <- layout(mat=matlayout, widths=c(1,.05, 1), 
         heights=c(.25,1, 1), FALSE)
par(plt=c(0,1,0,1));
plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1), log = "", asp = NA,
  xaxs="i", yaxs="i") 
  ##plot(x=c(0,1),y=c(0,1), type="l", main="", xlab="", ylab="",
  ##   xlim=c(.1,.9), ylim=c(.1,.9),
  ##       axes=FALSE);
  text(x=.15, y=.5, labels="Starting ", cex=2.5, pos=4);  
  text(x=.26, y=.5, labels=substitute( lambda[0] == a, 
       list(a=EM2LassoToy1b$LambdaAK[1]) ), cex=2.5, pos=4);
  text(x=.36,y=.5, labels=", ", cex=2.5, pos=4)
  text(x=.41,y=.5, labels=substitute( sigma == a,
        list(a=round(sqrt(EM2LassoToy1b$sigmaNoiseSq),3) )
      ),cex=2.5, pos=4
     )
   text(x=.53,y=.5, labels=", Multiplier = ", cex=2.5, pos=4)  
  text(x=.69,y=.5, labels=substitute( 3^a,
        list(a=
        round( log(EM2LassoToy1b$LambdaDK[2] / 
             EM2LassoToy1b$LambdaDK[1]) / log(3), 3 )  )
      ),cex=2.5, pos=4
     )
plot.new();
 plot.window(xlim=c(-.1,.1), ylim=c(0,1), log = "", asp = NA) 
  ##plot(x=c(0,1),y=c(0,1), type="l", main="", xlab="", ylab="",
  ##   xlim=c(.1,.9), ylim=c(.1,.9),
  ##       axes=FALSE);
  lines(x=c(0,0), y=c(0,1), lwd=3, col="black");
par(plt=c(.15,.95,.11,.85));
   headerWant = c(substitute( pi == a, list(a=EM2LassoToy1a$ppiuse)),""); 
   PrintSpreadFunction(EM2LassoToy1a, mfrows=FALSE,
      headerWANT = headerWant  );
      headerWant = 
        c(substitute( pi == a, list(a=EM2LassoToy1b$ppiuse)), "");
   PrintSpreadFunction(EM2LassoToy1b, mfrows=FALSE,
      headerWANT = headerWant  );  
      
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
