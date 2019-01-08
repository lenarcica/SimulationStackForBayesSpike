\name{Two2LassoCpp}
\alias{Two2LassoCpp}
\alias{Two2LassoCppFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Two2LassoCpp Function for producing Cpp objects as well as R5 Objects }
\description{
  Call Two2LassoCpp to do intense cross validation.  TwoLassoCpp is a simpler
  algorithm suitable for doing a single regression set. 
}
\usage{
Two2LassoCpp(X = -1, Y = -1, nSample = -1, 
  PiA = .5, SigmaSq = -999, 
  LambdaDK = -999, LambdaAK = -999, OrderSeq = -999,
  StLambdaD = -999, StLambdaA = -999, lambdaDMultC = 2^(2/8), 
  lambdaAMultC = 2^(-2/8), lambdaAmultStop = 20, TotalRuns = 20, 
  MaxCauchy = 20, CauchyEpsilon = .001, StartBeta = -999, StandardFlag = 0,
  XtX = -1, XtY = -1, InverseGammaConstant = 1, FixKa = -100, 
  InitKKs = -5, WLSWeights = -1, TDFNu = -1,
  Verbose = -1, SigmaVec = -1,
  RecordFlag = 0,
  PiAPrior = c(-1.5,-1.5), m1 = -1, m2 = -1, SigmaPrior = c(-1,-1),
  SigmaBar = -1, SigmaDf = -1, HoldOn = FALSE, Groupers = NULL,
  RunFlag = 1, L2ShrinkagePrior = c(-999,-999), 
  RecordL2Shrinkage = FALSE, DoCrossValidate = FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ Covariate matrix, dimensions n * p }
  \item{Y}{ Response vector, length n * 1 }
  \item{p}{ Number of Covariates}
  \item{XtX}{ Optional supply sum of squares covariates matrix }
  \item{XtY}{ Optional correlation of response to covariates vector}
  \item{dfTNoise}{ Noise level if to use TNoise }
  \item{n}{ Length of yys vector if not supplied }
  \item{pCoefs}{ Alias for p }
  \item{nSamples}{ Alias for n }
  \item{PiA}{ \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}}}{pi_A} probability  }
  \item{SigmaSq}{ \eqn{\sigma^2}{sigma^2} input value,  Note that if StandardFlag = 1
      Then USsigmaNoiseSq = sigmaNoiseSq / sd(yys) will be used in algorithm.}
  \item{PiAVector}{ \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}}}{pi_A} probability  vector if to use a sequence for all
    \eqn{\lambda_{\mbox{\tiny{$\mathcal{A}$}}}}{Lambda_A} values. }
  \item{SigmaSqVector}{ \eqn{\sigma^2}{sigma^2} input vector, if to use a sequence for all
    \eqn{\lambda_{\mbox{\tiny{$\mathcal{A}$}}}}{Lambda_A} values.}
  \item{PiAPrior}{A Beta distribution(2 paramter) Prior vector for
    \eqn{\pi_{\mbox{\tiny{$\mathcal{A}$}}}}{pi_A} probability  or a prior matrix
    if \code{LambdaAK} is along sequence. }
  \item{TotalRuns}{ Total number of iterations to increase \eqn{\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^{[\tau]}}{Lambda_D^[tau]}}
  \item{NumEMConv}{ Number of EM converge loops, usually 5-8 suffices }
  \item{MultEMCons}{ Allows one to reduce amount of EM loops if requested }
  \item{BetOlds}{ Input Beta vector if desired }
  \item{StandardFlag}{ Standardize Covariates and response if desired }
  \item{WLSWeights}{Vector of alternative weights for the data (in case of dfTNoise or Logit regression, set automatically) }
  \item{RunFlag}{ Whether to run algorithm after creating \code{\link{TwoLassoSEXP}} C++ object class }
  \item{DConfidenceInt}{ Not used anymore, configure Confidence interval estimation }
  \item{LambdaDK}{ Vector of \eqn{\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^[\tau]}{Lambda_D^[tau]} values if that is more practical }
  \item{LambdaAK}{  Vector of \eqn{\lambda_{\mbox{\tiny{$\mathcal{A}$}}}^[\tau]}{Lambda_A^[tau]} values if that is more practical }
  \item{OrderSeq}{ Sequence of integers length TotalRuns (+4,-4,5,-5,5...) For each value of 
     \eqn{\left(\lambda_{\mbox{\tiny{$\mathcal{D}$}}}^[\tau], 
        \lambda_{\mbox{\tiny{$\mathcal{A}$}}}^[\tau] \right)}{(Lambda_D^[t],Lambda_A^[t])} this decides whether to
         do E step first or M step (sign) and also how many EM steps per value (magnitude) }
  \item{NumCDOConv}{ Maximum number of times to run Coordinate Descent }
  \item{InitKKs}{ How much space to allocate intially for active coordinates. }
  \item{CDOEpsilon}{ Sum of deviations to quit Coordinate Descent Algorithms }
  \item{XTX}{ transpose of xx times xx if if desired }
  \item{XTY}{ transpose of xxs times yys if desired to pass less data }
  \item{InverseGammaConstant}{ Changes memory management within function when it inverts the gamma weights }
  \item{FixKa}{ If a positive nonzero value, implements Fermi-Dirac Approximation Limit Lasso, finding best step with FixKa active coefficients }
  \item{WLSWeights}{ If a positive vector of length yys, weights the data for Weighted Least Squares }
  \item{TDFNu}{ If a positive value, then ammends EM to weight Y data accoridng to t-distribution Noise of given coefficient.  Warning
   that in general, low degree of freedom T distributions produce multi-modal solutions, even in non-selection problems.  Limit-Lasso
   approach might be most helpful in generating large mode.}
  \item{RecordBetaCVFinish}{
    To do cross validation to solve for valid parameters, this part of TwolassoCpp must be allocated.
  }
  \item{DoCrossValidate}{Whether to do cross validation for this algorithm.}


}
\details{
   The TwoLassoCpp helper code algorithm for 2Lasso  Regression
}
\value{
    Returns \code{\link{TwoLassoSEXP}} object  
     \item{OnBeta}{Gives the final betas for all input values of LambdaA/LambdaD pairs}
     \item{BBOn1}{gives theestimates for expected B indicator values.}
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
 }
\author{ Alan Lenarcic }
\note{ 
    Limit-Lasso is obtained by giving LambdaAK/LambdaDK sequences,
     or alternatively, giving starting StLambdaA/StLambdaD and multipliers.
}
\seealso{ \code{\link{TwoLassoSEXP}}, \code{\link{EM2Lasso}}, \code{\link{CoordinateDescen}}, 
\code{\link{CoordinateDescentPrototype}} }
\examples{

library(TwoLassoCpp);
n = 200;
p = 10;
Beta = c(-1,1,-1,1,-1, rep(0,p-5));
X <- matrix(rnorm(n*p),n,p);
X <- StandardizeXX(X);

SigmaSq = .3;
Y <- X %*% Beta + sqrt(SigmaSq) * rnorm(n);
Y <- StandardizeXX(Y);
LambdaDK <- sqrt(n)*c(1,3, 4);
LambdaAK <- sqrt(n)* c(1,1/3,.25);
OrderSeq <- c(1,8,4);
PiA = .5;  CauchyEpsilon = .00001;  MaxCauchy = 200;
tt1 = NULL;  tt2 = NULL;  XtX = NULL;  XtY = NULL; StartBeta = rep(0,p);
dfTNoise = -1;  StandardFlag = 1;  ReturnBetas = NULL;  

AnsOb <- TwoLassoCpp(
     X = X, Y= Y, XtX = NULL, XtY = NULL,
     ReturnBetas = NULL, StartBeta = rep(0,p),
     dfTNoise = -1, StandardFlag = 1,
     CauchyEpsilon= .00001, MaxCauchy= 200,
     pCoefs = p, nSample = n, FixKa = -10,
     tt1 = NULL, tt2 = NULL, WLSWeights= NULL, 
     Verbose= 0, LambdaAK = LambdaAK, 
     LambdaDK= LambdaDK, OrderSeq= OrderSeq, PiA = .5, SigmaSq = SigmaSq,
     PiAVector= NULL, SigmaVector = NULL, PiAPrior = c(2,2),
     SigmaPrior= c(5,5),  InitKKs = 20, GroupSexp= NULL, 
     OnGammas = NULL, RecordFlag = 1, RecordBetaCVFinish = NULL,
     SigmaVectorInputs = NULL, PiAVectorInputs = NULL,
     LambdaAKVectorInputs= NULL, LambdADKVectorInputs = NULL,
     L2ShrinkagePrior = NULL, L2ShrinkageRecords = NULL,
     RunFlag = -1, RecordSigma = FALSE, TDFNu = 0, HoldOn = TRUE)
     
     
###### The Lars Analysis Format for the Paper
library(lars)
data(diabetes)
StlambdaA = 10; StlambdaD = StlambdaA;
XX <- StandardizeXX(diabetes$x);
YY <- StandardizeXX(diabetes$y);
MyLars <- lars(XX,YY);
## Member E is out
T1 <- proc.time();
TwoLassoToy <- TwoLassoCpp(X= XX, Y = YY, PiA = .5, 
  SigmaSq = .5, LambdaAK = c(1,.5,.1), LambdaDK = c(1,2,4),
  MaxCauchy = 7, MultEMCons = .99,  StandardFlag = 1,
  MaxCauchy = 200) 
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
          NumEMConv = 7, MultEMCons = .99, StandardFlag = 1) 
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
               NumEMConv = 7, MultEMCons = .99, StandardFlag = 1, 
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
               NumEMConv = 14, MultEMCons = .99, StandardFlag = 1, 
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
               NumEMConv = 14, MultEMCons = .99, StandardFlag = 1, 
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
