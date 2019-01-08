\name{0001AAAEfficientSimulator}
\alias{0001AAAEfficientSimulator.R}
\alias{0001AAAEfficientSimulator.r}
\alias{0001AAAEfficientSimulator}
\alias{EfficientSimulator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Efficient Simulator: similar defintiions for fixed effects estimators }
\description{
  Collection of Estimation Scripts for estimators considered for fixed effects sparse selection.
}
\usage{
library(TwoSimR5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\link{SMS}}{Input Simulation Object}
  
}
\details{
 
 File \code{0002AAB2LassoDefaultFunctions.r} largely gives the Function names used for all list of functions.  However a usual \code{\link{InputScriptFile}}
 should have a section redeclaring functions of interest for a given simulation set. However,
 please see \code{\link{0001AAAEfficientSimulator.r}} and \code{\link{0008AAAGroupEfficientSimulator.R}} for code that is detects elements of
 the simulation object "\code{\link{SMS}}" or result from "\code{\link{SimMeData}}()" and runs data extraction, and estimate extraction from the varied estimators.

 Note that hard criteria, like DIC, WLT, Cp, etc need to be developed for some estimators.  
 For 2Lasso and BayesSpike, much experimentation was developed over proper sparse priors.  A number of estimators were implemented first here (Horseshoe, Spike and Slab) 
 before a canoncial R-package version (Ishwaran Spike and Slab, monomvn BayesLassoa and Horseshoe) was discovered .  While we tried our best, for the most part the CRAN version of 
 these estimators is superior, and a better test of an estimator's performance.
 
 
 Note that many more estiators exist in this file than were featured in papers.  This is because some additional settings and estimators were deemed redundant, 
 that some estimation methods performed so poorly as to not be competitive.  It is likely that mistakes were made during our implementation of each method, but the best we can do is test to demonstate a consistent perfomrance, share the code used to in the procedure, and hope that
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
SMS <- SimMeData();
MCPEstimate <- GenerateMCP(SMS, LARSSeekFlag, DoConvexMin=1, DoCI = DefaultDoCI)
## q();
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ts }% __ONLY ONE__ keyword per line
