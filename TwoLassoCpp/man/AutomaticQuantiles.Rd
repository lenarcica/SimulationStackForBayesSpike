\name{AutomaticQuantiles}
\alias{AutomaticQuantiles}
\alias{automaticquantiles}
\alias{ConfidenceMatrix}
\alias{ConfidenceQuantiles}
\alias{ConfidenceQuartiles}
\alias{CredibilityMatrix}
\alias{Credibility}
\alias{credibility}
\alias{confidencequantiles}
\alias{confidencematrix}
\alias{confidencequartiles}
\alias{credibilityquantiles}
\alias{CredibilityQuantiles}
\alias{CredibilityMatrix}
\alias{quantiles}
\alias{Quantiles}
\alias{automaticQuantiles}
\alias{SymmetricQuantiles}
\alias{HPDQuantiles}
\alias{BetaSymmetricQuantiles}
\alias{BetaSymmetricUnshrinkQuantiles}
\alias{BetaSymmetricUnshrunkQuantiles}
\alias{BetasymmetricQuantiles}
\alias{betaSymmetricQuantiles}
\alias{BetaHPDUnshrinkQuantiles}
\alias{BetaQuantiles}
\alias{betaQuantiles}
\alias{QuantilesBeta}
\alias{BetaHPDQuantiles}
\alias{betaHPDQuantiles}
\alias{Betahpdquantiles}
\alias{BetaHpdQuantiles}
\alias{betaHpdQuantiles}
\alias{TauSymmetricQuantiles}
\alias{TausymmetricQuantiles}
\alias{tauSymmetricQuantiles}
\alias{TauQuantiles}
\alias{tauQuantiles}
\alias{QuantilesTau}
\alias{TauHPDQuantiles}
\alias{tauHPDQuantiles}
\alias{Tauhpdquantiles}
\alias{TauHpdQuantiles}
\alias{tauHpdQuantiles}
\alias{PiASymmetricQuantiles}
\alias{piasymmetricQuantiles}
\alias{piASymmetricQuantiles}
\alias{PiAQuantiles}
\alias{piAQuantiles}
\alias{QuantilesPiA}
\alias{PiAHPDQuantiles}
\alias{piaHPDQuantiles}
\alias{PiAhpdquantiles}
\alias{PiAHpdQuantiles}
\alias{piAHpdQuantiles}
\alias{SigmaSymmetricQuantiles}
\alias{SigmasymmetricQuantiles}
\alias{sigmaSymmetricQuantiles}
\alias{sigmaQuantiles}
\alias{sigmaQuantiles}
\alias{QuantilesSigma}
\alias{SigmaHPDQuantiles}
\alias{sigmaHPDQuantiles}
\alias{Sigmahpdquantiles}
\alias{SigmaHpdQuantiles}
\alias{SigmaHpdQuantiles}
\alias{SigSymmetricQuantiles}
\alias{sigSymmetricQuantiles}
\alias{SigHPDQuantiles}
\alias{SighpdQuantiles}
\alias{KeepPosteriorQuantiles}
\alias{keepPosteriorQuantiles}
\alias{KeepQuantiles}
\alias{keepquantiles}
\alias{keepposteriorquantiles}
\alias{Keepposteriorquantiles}
\alias{KeepposteriorQuantiles}
\alias{KeepPostQuantiles}
\alias{keepPostQuantiles}
\alias{LambdaIndexForConfidenceIntervals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Marginal Posterior Quantiles of Draws}
\description{

   Generte Quantiles of the Marginal posterior of 2 Lasso distribution. 
   
   To calculate, one must set "\code{ConfidenceQuantiles}" element of 
     \code{\link{TwoLassoSEXP}} C++ object class.  Then, the element
     \code{ConfidenceMatrix} can be taken from TwoLassoSEXP assuming the algorithm
     has been run past the point set for \code{LambdaIndexForConfidenceIntervals}
   
   An example of quantiles to look for consider:
   
   \code{TwoLassoSEXPOb$ConfidenceQuantiles <- c(.5, .25,.75, .05,95,.025,.975,.01,.99);}
   
   
   
  
   
}
\usage{
TwoLassoSEXPOb$ConfidenceQuantiles <- c(.5, .25,.75, .05,95,.025,.975,.01,.99);
TwoLassoSEXPOb$LambdaIndexForConfidenceIntervals <- 1;
TwoLassoSEXPOb$ConfidenceMatrix; 
TwoLassoSEXPOb$UnshrinkConfidenceMatrix; 
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{ConfidenceQuantiles}{Set a vector of quantiles in the of \code{\link{TwoLassoSEXP}} object.  
    Set to determine what intervals are recorded.  Numbers should be within \code{[0,1]} }
  \item{LambdaIndexForConfidenceIntervals}{An index of the LambdaD sequence to find point estimate of beta to use to calculate object}


}
\details{
   
}
\value{

  \item{\code{ConfidenceMatrix}}{ Output matrix that delivers quantiles of marginal posterior given LambdaA and LambdaD}
  \item{\code{UnshrinkConfidenceMatrix}}{ Removes selection from Credibility interval for parameters.  Should have
    better Confidence frequency performance near zero.  }
 
}
\references{ TwoLasso Thesis }
\author{ Alan Lenarcic}
\note{ 
}
\seealso{ \code{\link{TwoLassoRegression}}, \code{\link{TwoLassoSEXP}}
  }
\examples{
library(TwoLassoCpp);
n = 200;
p = 10;
Beta = c(-1,1,-1,1,-1, rep(0,p-5));
X <- matrix(rnorm(n*p),n,p);
X <- StandardizeXX(X);

SigmaSq = .3;
Y <- X %*% Beta + sqrt(SigmaSq) * rnorm(n);
LambdaDK <- sqrt(n)*c(1,3, 4);
LambdaAK <- sqrt(n)* c(1,1/3,.25);
OrderSeq <- c(1,8,4);
PiA = .5;  CauchyEpsilon = .00001;  MaxCauchy = 200;
tt1 = NULL;  tt2 = NULL;  XtX = NULL;  XtY = NULL; StartBeta = rep(0,p);
dfTNoise = -1;  StandardFlag = 1;  ReturnBetas = NULL;  

AnsOb <- TwoLassoCpp(
     X = X, Y= Y, XtX = NULL, XtY = NULL,
     ReturnBetas = NULL, StartBeta = rep(0,p),
     dfTNoise = -1, StandardFlag = 0,
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
AnsOb$TLS$ConfidenceQuantiles <- c(.5, .25, .75, .05,.95, .025, .975, .01, .99, .005, .995);
AnsOb$TLS$Verbose = 5; 
AnsOb$TLS$ConfidenceMatrix;
AnsOb$TLS$UnshrinkConfidenceMatrix;

TC <- 100;
CovInt = matrix(0, TC, 5);
WInt = matrix(0, TC, 5);
OnCovInt = matrix(0, TC, 5);
OnWInt =  matrix(0, TC, 5);
OffCovInt = matrix(0, TC, 5);
OffWInt =  matrix(0, TC, 5);

n = 50;
p = 75;

ii <- matrix(rep(1:p, each=p), p, p);
jj <- t(ii);
corrM <- .2^(abs(ii - jj))
LCorr <- t(chol(corrM));
MyBad <- matrix()
for (ii in 1:TC) {

Beta = .2*c(-1,1,-1,1,-1,1, rep(0,p-6));
X <- matrix(rnorm(n*p),n,p) %*% LCorr;
X <- StandardizeXX(X);

SigmaSq = .5;
Y <- X %*% Beta + sqrt(SigmaSq) * rnorm(n);
LambdaDK <- sqrt(n)*c(1,3, 4, 8);
LambdaAK <- sqrt(n)* c(1,1/3,.25, .05);
OrderSeq <- c(1,8,4,4);
PiA = .5;  CauchyEpsilon = .00001;  MaxCauchy = 200;
tt1 = NULL;  tt2 = NULL;  XtX = NULL;  XtY = NULL; StartBeta = rep(0,p);
dfTNoise = -1;  StandardFlag = 1;  ReturnBetas = NULL;  

AnsOb <- TwoLassoCpp(
     X = X, Y= Y, XtX = NULL, XtY = NULL,
     ReturnBetas = NULL, StartBeta = rep(0,p),
     dfTNoise = -1, StandardFlag = 0,
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
ABStart <- c(.5, .25, .75, .05,.95, .025, .975, .01, .99, .005, .995);
ABCorrect <- (pmin(ABStart, 1-ABStart) * (n/(n+log(p))))
ABCorrect[ABStart > .5] = 1-ABCorrect[ABStart>.5];
ABCorrect[ABStart == .5] = .5;
ABCorrect = ABStart;
AnsOb$TLS$ConfidenceQuantiles <- ABCorrect;
AnsOb$TLS$LambdaIndexForConfidenceIntervals <- 3;
AnsOb$TLS$Verbose = 5; 
MyMs <- t(matrix(2:11,2, 5));
CCs <- c(50, .9, .95, .98, .99);
AC <- AnsOb$TLS$ConfidenceMatrix;
for (jj in 1:length(CCs)) {
  WInt[ii,jj] = mean(abs(AC[,MyMs[jj,2]] - AC[,MyMs[jj]]));
  MyCo <- rep(0, p);
  MyCo[AC[,MyMs[jj,2]] >= Beta & AC[,MyMs[jj,1]] <= Beta] <- 1;
  CovInt[ii,jj] <- mean(MyCo);
  OnWInt[ii,jj] = mean(abs(AC[Beta != 0.0 ,MyMs[jj,2]] - AC[Beta != 0.0,MyMs[jj]]));
  MyCo <- rep(0,length(Beta[Beta != 0.0]) );
  MyCo[AC[Beta != 0.0,MyMs[jj,2]] >= Beta[Beta != 0.0] & AC[Beta != 0.0 ,MyMs[jj,1]] <= Beta[Beta != 0.0]] <- 1;
  OnCovInt[ii,jj] <- mean(MyCo);
  OffWInt[ii,jj] = mean(abs(AC[Beta == 0.0,MyMs[jj,2]] - AC[Beta == 0.0,MyMs[jj]]));
  MyCo <- rep(0, length(Beta[Beta==0.0]) );
  MyCo[AC[Beta == 0.0,MyMs[jj,2]] >= 0.0 & AC[Beta == 0.0,MyMs[jj,1]] <=  0.0] <- 1;
  OffCovInt[ii,jj] <- mean(MyCo);
  
}
}

colMeans(OffCovInt);  colMeans(OnCovInt);  colMeans(CovInt);
colMeans(WInt); colMeans(OnWInt);  colMeans(OffWInt);


piA = .4;  A = 4; B = 2; LambdaA = 1; LambdaD = 4;
AnsOb$TLS$TestAllInt(A,B,piA,LambdaA,LambdaD);
AnsOb$TLS$TestSuperInt(A,B,piA,LambdaA,LambdaD, 
  AnsOb$TLS$TestAllInt(A,B,piA,LambdaA,LambdaD), 0)
XX <- (-500:500)/500 * 10;
DRX <- rep(0, length(XX));
ARX <- rep(0, length(XX));
DRX <- exp(A * XX - B * XX^2/2) * ( piA * exp(-LambdaA * abs(XX)) + (1-piA) * exp(-LambdaD * abs(XX)));
for (ii in 1:length(XX)) {
  ARX[ii] <- AnsOb$TLS$TestSuperInt(A,B,piA,LambdaA,LambdaD, 
    AnsOb$TLS$TestAllInt(A,B,piA,LambdaA,LambdaD), XX[ii]);
}
par(mfrow=c(1,2));
plot(DRX~XX, type="l");  plot(ARX~XX, type="l");

ARD <- AnsOb$TLS$TestSeekQuantile(.99999, -1, .0000001,
    A, B, PiA, LambdaA, LambdaD, AnsOb$TLS$TestAllInt(A,B,PiA,LambdaA, LambdaD))
    
AnsOb$TLS$TestSuperInt(A,B,piA,LambdaA,LambdaD, 
  AnsOb$TLS$TestAllInt(A,B,piA,LambdaA,LambdaD), ARD);

     
Two2LassoCpp <- function(X = X, Y = Y, nSample = n, 
  PiA = PiA, SigmaSq = SigmaSq, 
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ts }% __ONLY ONE__ keyword per line
