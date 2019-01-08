\name{RunMakeScriptFiles}
\alias{RunMakeScriptFiles}
\alias{RunMakeScriptFiles.r}
\alias{runmakescriptfiles}
\alias{Runmakescriptfiles}
\alias{ScriptGetArgs}
\alias{CheckUniqueProcessIdentifier}
\alias{MakeDirSaves}
\alias{UpdateNameFunctions}
\alias{AssessFirstTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ RunMakeScriptFiles }
\description{
  The process of generating new simulations based upon a \code{\link{InputScriptFile}}, identifying a save directory, saving the
  simulation data, input parameters in a \code{\link{SimMeData}} object, and assessing how many completed simulation estimates exist in 
  the database.  When every estimator has fit every simulation (say 1000 simulations each for approx 55 different estimators),
  this script will stop contributing more simulations or estimates.  After this point \code{\link{CollectiSimulationResults}}
  is ready to run. 
}
\usage{
library(TwoSimR5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\link{InputScriptFile}}{A model input file communicating size, p, etc choices}
  \item{\link{RunMakeScriptFiles}}{The step of generating simulations in parallel, and assessing recorded simulation results database}
  \item{ScriptGetArgs}{Code that is run after the InputScriptFile to determine 
    all remaining undeclared arguments.  Note that a complete InputScriptFile potentially covers much information, 
    including which statistical estimators to use and settings to choose}
  \item{AssessFirstTime}{This code actually assess how much simulations still have to be simulated, 
  and how many estimators need to be tested on each simulation.  Unfortunately for it's bad name, 
  it is run everytime.  It is a part of \code{2LassoMaster} file.}
}
\details{
 Some attempt was made to use "R5" objects to improve object oriented coordination in this module.
 
 One of the mose important objects is "\code{\link{SMS}}" or result from \code{\link{SimMeData}}() which stores a simulations object and relevant
 details which all of the simulation estimators collectively take in.
 
 Much data is saved to the global environment to generate these simulations.
}
\value{
    
    Major Analysis is usually a single vector Beta (usually the model after fitting a penalty-free version on the active set
    determined by the algorithm, or a Median model from Gibbs Samplers), Model Inclusion Probabilities if given,
    identies of rhe active-chosen variables, results of quantiles of credibility intervals coverage, losses like l2 loss, sphere loss
    on the fit.
    
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
## A Script demonstrating steps to  RunMakeScriptFiles.r: 
## R --vanilla --args PracticeTestn100p100.r 1 1 1 1 1 < RunMakeScriptFiles.r
## R --vanilla --args GroupLongTestn100p10000.r 1 1 1 1 1 
library(TwoSimR5);    library(AlanDirectories);

TwoSimR5:::ScriptGetArgs();
CheckUniqueProcessIdentifier();    ## In 2LassoSaveRoutines.$

try(library(TwoLassoCpp, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);

verbose = 1;
Practice = 1;
if (length(list.files("c:/Stat")) >= 1) {
  argeese <- paste(DropHome, "//TwoLassoPackage//PracticeTestn100p100.r", sep="");
}
source(paste(as.character(argeese), sep=""));
                   

TwoSimR5:::MakeDirSaves();
JustPrintFile <- "MakingASimulation.txt";
UpdateNameFunctions();

NowDelete=TRUE;
AssessFirstTime(ALargeContainDir, ASmallContainDir, TargetTotal = TargetTotal, TotalThralls=TotalThralls);


##bash ~/TR5Script//Test2013April//N100p1000k6sigp5piSV3sigSV3/RBayesAlgorithm1.txt


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ts }% __ONLY ONE__ keyword per line
