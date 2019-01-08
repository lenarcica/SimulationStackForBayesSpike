\name{GetForMeAll2}
\alias{GetForMeAll2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetForMeAll2(ppiuseV, sigmaFactV, ppCovV, NNOnV, NROnV, NoiseOnV, NumReps, PrintFlag = 0, BRealInput = c(1, -1, 1, -1, 1, -1), LARSSeekFlag = 0, musize = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppiuseV}{ Unused Active variable}
  \item{sigmaFactV}{ Mean and variance of active correlations between covariates }
  \item{ppCovV}{ Probability of active correlation between two covariates }
  \item{NNOnV}{ Number of Samples per simulation }
  \item{NROnV}{ Total Number of Covariates }
  \item{NoiseOnV}{ Noise level of data }
  \item{NumReps}{ Number of repetitions (200 default) }
  \item{PrintFlag}{ Whether to print information while simulating }
  \item{BRealInput}{An Input vector to start seek for the active variables  if desired.}
  \item{LARSSeekFlag}{ A factor to begin Two-Lasso with relative to Lasso hit. }
  \item{musize}{ Length of true active vector. }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  In one line of OneVV input first 
     #1  Largely anachronistic integer
     #2  Deviation of correlation between covariates 
     #3  Probability of correlation between covariates
     #4  Total # of N samples
     #5  Total # of K covariates
     #6  Sigma Noise Level
}
\references{ Alan Lenarcic's 2009 Harvard Thesis }
\author{ Alan Lenarcic }
\note{ Enjoy
}
\seealso{ GetForMeAAA2 }
\examples{

 ##################################################################3
  ####################   Table 1   
  ####################DefinitivePrint;   
  #### Second Paper ##6
  PrMeVec <- c(4,7,9,8,10); 
  ##GFMAA <- TRAA2;
  ##OneVV <- OneVV;
  NCount <- 200; KPAm <- 6;
   TMM = FALSE; rndit = 2;
 ### Noise1
   OneVV <-  rbind(  
                       c(7, .2, .3, 25, 20, 1),
                       c(7, .2, .3, 25, 30, 1),
                       c(7, .2, .3, 25, 40, 1),
                       c(7, .2, .3, 10, 25, 1),
                       c(7, .2, .3, 15, 50, 1),
                       c(7, .2, .3, 25, 50, 1),
                       c(7, .2, .3, 35, 50, 1),
                       c(7, .2, .3, 45, 50, 1),
                       c(7, .2, .3, 25, 100, 1),
                       c(7, .2, .3, 50, 100, 1),
                       c(7, .2, .3, 75, 100, 1),
                       c(7, .2, .3, 100, 100, 1));  
 ## GFMAA <-   GetForMeAAA2(  OneVV, NumReps = 200, PrintFlag = 0, BRealInput = c(1,-1,1,-1,1,-1),  LARSSeekFlag = 2.5 )                                          
  

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datasets }
\keyword{ methods }% __ONLY ONE__ keyword per line
