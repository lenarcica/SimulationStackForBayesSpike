\name{SimMeData}
\alias{SimMeData}
\alias{simmedata}
\alias{SimMeData.R}
\alias{SimGroupData}
\alias{simgroupdata}
\alias{SMS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate Simulation Objects }
\description{
This file takes a number of possible covariance structures and generates data from a linear model, along with necessary n,p,etc characteristics
useful for conducting a standardized structure necessary for testing sparse fixed effects and sparse grouped effects models.
}
\usage{
SimMeData(n = -999, p = -999, k = -999, sigma = -999,
     NLen = -999, kLen = -999, kActiveLen = -999,  ElseGetX = NULL, 
     SigmaNoise = -999, GenerateBetaVec = DefaultGenerateBetaVec, 
     CorrelationXmatrix = DefaultCorrelationXmatrix, tNoiseDF = DefaultNoiseDF,
     LogitRegression = FALSE, LogitFactor = 1, Beta0 = 0, ExperimentName="", jobii= 0, 
     WorkingRow = 0, AlreadyLocked=TRUE, MeanCenterXX = TRUE,
     UniqueProcessIdentifier = "", ISample = 1,
     DontSave=FALSE, InOrder=FALSE, dfRobit = -1)
     
SimGroupData(n = -999, pGroups = -999, GroupSize=5, kGroups = -999, sigma = -999,
     SigmaNoise = -999, tNoiseDF = -1, GenerateBetaGroupVec = DefaultGenerateBetaGroupVec, 
     CorrelationXmatrix = DefaultCorrelationXmatrix, ElseGetX = NULL,  ElseTauEndList = NULL,
     LogitRegression = FALSE, LogitFactor = 1, Beta0 = 0, ExperimentName="", jobii= 0, 
     WorkingRow = 0, AlreadyLocked=TRUE, MeanCenterXX = TRUE,
     UniqueProcessIdentifier = "", ISample = 1, GroupsSumToZero=0,
     DontSave=FALSE, InOrder=FALSE, dfRobit=-1)
     
DoubleSimulation(SMS, AlreadyLocked, DontSave=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{count independent datapoints}
  \item{p}{count covariates}
  \item{k}{total count of non zero covariages}
  \item{sigma}{standard deviation of noise, will be squared into SigmaSqNoise}
  \item{GenerateBetaVec}{A function describing how to generate nonzero Betas of length k}
  \item{CorrelationXmatrix}{A function describing how to get correlation structure of interest}
  \item{tNoiseDF}{A number of degrees of freedom for tNoise}
  \item{LogitRegression}{Generate a Logit bernoulli result for each covariate?}
  \item{LogitFactor}{If generating Logistic regression, spread data by a given factor}
  \item{Beta0}{Logit regression, what should be the intercept?}
  \item{MeanCenterXX}{Remove Mean and center from X covariates?}
  \item{UniqueProcessIdentifier}{What will be assigned as a unique process that generated this data?}
  \item{DontSave}{Save this SMS object to relevant saves directory or not?}
}
\details{
Will generate a "Noised" and "UnNoised" copy of itself.  The Noised version does not contain accurate k or sigma information.  Some estimators will be 
given both a noised and unnoised version of this material to see how a vague amount of prior information on sparsity informs on performance.

\code{SimGroupData()} simulates grouped Betas for grouped data.  Generating a covariance structure for correlated groups is a little awkward, so we did not set a hard rule on forcing this for our grouped simulation (this would have to make 
nearby groups correlated on one or more group factors, and is a bit too open, or genetics-specific a question to model for a statistics paper.)

\code{DoubleSimulation}() simply takes an existing SMS object and creates a copy of existing X columns in the next p columns of an object.  This code is only necessary for testing perfect repeated modes in the data.
}
\value{
    
    This SMS is a submitted simulation object for which simulation analysis can be performed.
    
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
