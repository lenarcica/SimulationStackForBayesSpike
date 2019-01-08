\name{0008AAAGroupEfficientSimulator}
\alias{0008AAAGroupEfficientSimulator.R}
\alias{0008AAAGroupEfficientSimulator.r}
\alias{0008AAAGroupEfficientSimulator}
\alias{GroupEfficientSimulator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Group Efficient Simulator: similar defintiions for fixed effects estimators }
\description{
  Collection of Estimation Scripts for estimators considered for fixed effects sparse selection.   Grouped versions from \link{EfficientSimulator} are featured in this file
}
\usage{
library(TwoSimR5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\link{SMS}}{Input Simulation Object}
  \item{pGroups}{Number of Groups instead of number of paramters}
  
}
\details{
 
 File \code{0002AAB2LassoDefaultFunctions.r} largely gives the Function names used for all list of functions.  However a usual \code{\link{InputScriptFile}}
 should have a section redeclaring functions of interest for a given simulation set. However,
 please see \code{\link{0001AAAEfficientSimulator.r}} and \code{\link{0008AAAGroupEfficientSimulator.R}} for code that is detects elements of
 the simulation object "\code{\link{SMS}}" or result from "\code{\link{SimMeData}}()" and runs data extraction, and estimate extraction from the varied estimators.

 Note that grouped data simulation presents different challenges from un grouped, and defining a realistic correlation structure amongst nearby groups was deemed too 
  domain specific.  

 For 2Lasso and BayesSpike, much experimentation was developed over proper sparse priors. 
 
}
\value{
    
   Each function listed here returns a \code{SubMitList} Object
   Statistics collected are Beta L2, Sphere loss, Type1, Type2, timing, etc + MIP performance and Credibility Interval performance if applicable.
   
   Example Return Structure:
   \item{type}{Example ="TwoLassoCV", name of estimator used}
	 \item{BetaFit}{Best Beta vector of length p returned, to be used for L2 Error}
   \item{BBHatsFit}{Indicator which elements of Beta are really chosen as active}
   \item{FitTime}{Fit time of the algorithm, minus some busywork}
	 \item{OtherBetas}{If other Betas, Mean Posterior instead of Median, fit under penalty as opposed to fit of active after penalty to consider}
   \item{EMLARSObject}{Object not usually supplied}
   \item{CIEst}{Confidence Interval locations}
   \item{CITime}{Time to estimate Confidence Intervals}
   \item{CIQuantiles}{Confidence Interval Quantiles .50, .75, .95 etc... estimated}
   \item{MIPReport}{Model Inclusion Probabilities if reported}
   \item{DoLogit}{Is this an estimate of a Logit Model?}
       
}
\references{  BayesSpike submission, Perhaps a TwoLasso Whitepaper }
\author{ Alan Lenarcic, Will Valdar }
\note{ 

}
\seealso{ 
  \code{\link[BayesSpike]{BayesSpikeRegression}}, 
   \code{\link[code]{mcmc}},
  \code{\link[methods]{ReferenceClasses}},
  \code{\link{Save}}, \code{\link{ReLoad}},
  }
\examples{
## An example Running the CollectSimulationsResults process:
## Opening to RunMakeScriptFiles.r
## R 
## R --vanilla --args LongTestn100p1000.r 1 1 1 1 1 1 1 1 
library(TwoSimR5);    library(AlanDirectories);
n <- 100; p <- 25; sigma <- 1;
SMS <- SimGroupData();
GRLEstimate <- GenerateGroupLassoReg(SMS, Interval = Interval,
  Criterion = "BIC")
## q();
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ts }% __ONLY ONE__ keyword per line
