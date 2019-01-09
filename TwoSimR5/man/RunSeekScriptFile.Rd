\name{RunSeekScriptFile}
\alias{RunSeekScriptFile}
\alias{RunSeekScriptFile.r}
\alias{RunAProcessExperiment}
\alias{Runaprocessexperiment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ RunAProcessExperiment }
\description{
  Onces simulations have been generated, this code exectues an estimate from the missing, uncompleted simulations and algos.  Although a single algo could fail on everythese, ideally most estimates will complete and the problematic estimator can be identified and improved.
  algos will result and can be recorded 
}
\usage{
library(TwoSimR5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\link{InputScriptFile}}{A model input file communicating size, p, etc choices same InputLinkScriptFile is used to determine settings.  Note this should be the same file pased to a RunMakeScriptFiles.r call first both befoe and at end.}
  \item{\link{RunSeekScriptFiles}}{The step of running estimatlrs for estimators in parallel, and assessing recorded simulation results database}
  \item{ScriptGetArgs}{Code that is run after the InputScriptFile to determine 
    all remaining undeclared arguments.  Note that a complete InputScriptFile potentially 
    covers much information, including which statistical estimators to use and which choices to make}
  \item{RunAProcessExperiment}{This code runs estimator run and conducts basic testss. It is a part of \code{2LassoMaster}.}
}
\details{
 The "RunSeekScriptFile" is an R script that prepares an R instance opened as a thread in LSF.  This thread will identify which
 Estimator and simulation it is assigned.  It then runs \code{RunAProcessFile()} on this combo to generate a fit and then asses the performance.
 
 
 Some attempt was made to use "R5" objects to improve object oriented coordination in this module.
 
 Please see \code{\link{0001AAAEfficientSimulator.r}} and \code{\link{0008AAAGroupEfficientSimulator.R}} for code that is detects elements of
 the simulation object "\code{\link{SMS}}" or result from "\code{\link{SimMeData}}()" and runs data extraction, and estimate extraction from the varied estimators.

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
## An example Running the RunSeekScriptFile process:
##R --vanilla --args LongTestn100p1000.r 1 1 1 1 1 1 1 1 
library(TwoSimR5);

## Preamble detects which simulation and estimator we should work on next.
eval(parse(text=RunSeekScriptPreamble()));
if (FALSE) {
  ## How to replace RunSeekScriptPreamble with a search for a particular estimator of interest.
  ABB <- (1:NROW(SubMitList))[SubMitList[,3] == "0" & SubMitList[,2]=="GenerateBayesSpikeAutoPrior"][1]
  SMSName=SubMitList[ABB,1];
  FunctionName =SubMitList[ABB,2];
  DoCI = DefaultDoCI;
}

## Run the analysis.
ABack <- TwoSimR5:::RunAProcessExperiment(SMSName=SubMitList[ABB,1], 
  FunctionName =SubMitList[ABB,2], verbose = verbose, DoCI = DefaultDoCI)
SubMitList[ABB,3] <- as.character(ABack);
ABFinish <- proc.time();
if (ABFinish[3] - ABStart[3] <= 5) {
  Sys.sleep(5- ABFinish[3]+ABStart[3])
}
setwd(ASmallContainDir);

## Save result in submissions directory.  If SubMitList has confidence/credibility interval
## estimates, it will need to save this larger information in ALargeContainDir.
try(save(SubMitList=SubMitList, ASmallContainDir=ASmallContainDir, 
       ALargeContainDir=ALargeContainDir, 
       CurrentLargeContainDir=CurrentLargeContainDir, 
       CurrentSmallContainDir=CurrentSmallContainDir, file=paste("Submissions//", WantFile, sep="")));
## q();
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ts }% __ONLY ONE__ keyword per line
