\name{CollectSimulationResults}
\alias{CollectSimulationResults.R}
\alias{collectsimulationresults}
\alias{CollapseOnFunction}
\alias{collapseonfunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Collapse simulations into one RData table }
\description{
  Onces simulations and estimators are generated, collect the temporary RData saves into one file.
}
\usage{
library(TwoSimR5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\link{OnFunctionNum}}{Should Detect the list of DefaultUseFunctions used in this simulation.  Open all of the save files
   for this simulation and combine into a single table useful for taking performance statistics}
  
}
\details{
 Some attempt was made to use "R5" objects to improve object oriented coordination in this module.
 
 File \code{0002AAB2LassoDefaultFunctions.r} largely gives the Function names used for all list of functions.  However a usual \code{\link{InputScriptFile}}
 should have a section redeclaring functions of interest for a given simulation set. However,
 please see \code{\link{0001AAAEfficientSimulator.r}} and \code{\link{0008AAAGroupEfficientSimulator.R}} for code that is detects elements of
 the simulation object "\code{\link{SMS}}" or result from "\code{\link{SimMeData}}()" and runs data extraction, and estimate extraction from the varied estimators.

}
\value{
    
   Statistics collected are Beta L2, Sphere loss, Type1, Type2, timing, etc + MIP performance and Credibility Interval performance if applicable.
    
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

TwoSimR5:::ScriptGetArgs();
CheckUniqueProcessIdentifier();    ## In 2LassoSaveRoutines.$

verbose = 1;
Practice = 1;
source(paste(as.character(argeese), sep=""));

CheckArgFunction();                    

TwoSimR5:::MakeDirSaves();

UpdateNameFunctions();
JustPrintFile <- paste("CollectingASimulation", LTO, ".txt", sep="");

for (iOni in OnFunctionNums) {
  ## Collapse for every OnFunction in all functions used for this simulation
  try(CollapseOnFunction(OnFunction = NameFunctions[iOni]));  
}
## q();
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ts }% __ONLY ONE__ keyword per line
