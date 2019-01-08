\name{TwoSimR5}
\alias{TwoSim}
\alias{TwoSimR5}
\alias{TwoSimR}
\alias{TwoSimr} 
\alias{twosimr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ TwoSimR5 }
\description{
  A Simulation settings package designed around many linear, grouped, sparse estimators.
  
  Most R packages for sparse selection unfortunately have different requirements and outputs.
  
  The purpose of this package was to aid in a simulation set generating simulations from a common
  framework (including simulation of noised prior information), recording the choices of that
  framework, assessing the need for remaining simulations, assessing simulation failure, and
  generating a series of tests on each parameter model fit.
  
  At this point, the vast majority of estimators come from 3rd party packages considered
  canoncial versions of a simulation method.  But many selection criteria
  (for instance WLT for SCAD), or DIC, needed to generated in code.
}
\usage{
library(TwoSimR5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\link{InputScriptFile}}{A model input file communicating size, p, etc choices}
  \item{\link{RunMakeScriptFiles}}{The step of generating simulations in parallel, and assessing recorded simulation results database}
  \item{\link{RunSeekScriptFile}}{Run a simulation fit from one estimator abd save result}
  \item{\link{CollectSimulationResults}}{Colleect the separate simulation estimation results files into common matrices for averages/sds}
  \item{\link{0001AAAEfficientSimulator}}{This file contains all setting sand coding choices for fixed effects sparse estimators}
  \item{\link{0008AAAGroupEfficientSimulator}}{This file contains all settings and coding choices for "Grouped Estimators"}
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
try(require(TwoSimR5, quietly=TRUE);
## We should give an example of full automation process.
## However, users would need to taylor simulation process to their own server
## system.
## TwoSimR5 has been developed with a "lsf" system in mind, though there are alternatives.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ts }% __ONLY ONE__ keyword per line
