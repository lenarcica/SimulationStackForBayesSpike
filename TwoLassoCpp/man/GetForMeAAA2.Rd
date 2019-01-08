\name{GetForMeAAA2}
\alias{GetForMeAAA2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ GetForMe Simulations }
\description{
  Functions for simulating data and estimator 
        fits and rating success for the data
}
\usage{
GetForMeAAA2(OneVV, NumReps, PrintFlag = 0, BRealInput = c(1, -1, 1, -1, 1, -1), LARSSeekFlag = 0, musize = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OneVV}{ Matrix of Parameter Input values, there are 5 inputs per sim. }
  \item{NumReps}{ NumReps is total number of simulations to save, usually 200 }
  \item{PrintFlag}{ If positive integer, prints info to console about progress }
  \item{BRealInput}{ Real Values for input variable }
  \item{LARSSeekFlag}{ Describes consant to divide LARS-fixed lambda by to use as Limit-Lasso Input }
  \item{musize}{ The Length of BRealInput }
}
\details{
  In one line of OneVV input first 
     #1  Largely anachronistic integer
     #2  Deviation of correlation between covariates 
     #3  Probability of correlation between covariates
     #4  Total # of N samples
     #5  Total # of K covariates
     #6  Sigma Noise Level
}
\value{
    Returns Sim List
}
\references{ R LARS package, Efron  and Tibshirani.  Lin and Yuan: "Lasso w=1"}
\author{ Alan Lenarcic }
\note{   Enjoy
}
\seealso{ GetForMeAAA }
\examples{
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
  ##GFMAA <-   GetForMeAAA2(  OneVV, NumReps = 200, PrintFlag = 0, BRealInput = c(1,-1,1,-1,1,-1),  LARSSeekFlag = 2.5 )                                          
  ##A1 <- DoAllTheSaving(GFMAA, OneVV, KPAm, PrMeVec, NCount)
  ##A1b <- DoAllTheSavingL(GFMAA, OneVV, KPAm, PrMeVec, NCount)
   OneVV <-  rbind(  
                       c(7, .2, .3, 25, 20, .5),
                       c(7, .2, .3, 25, 30, .5),
                       c(7, .2, .3, 25, 40, .5),
                       c(7, .2, .3, 10, 25, .5),
                       c(7, .2, .3, 15, 50, .5),
                       c(7, .2, .3, 25, 50, .5),
                       c(7, .2, .3, 35, 50, .5),
                       c(7, .2, .3, 45, 50, .5),
                       c(7, .2, .3, 25, 100, .5),
                       c(7, .2, .3, 50, 100, .5),
                       c(7, .2, .3, 75, 100, .5),
                       c(7, .2, .3, 100, 100, .5));  
  ##GFMAA <-   GetForMeAAA2(  OneVV, NumReps = 200, PrintFlag = 0, BRealInput = c(1,-1,1,-1,1,-1),  LARSSeekFlag = 2.5 )                                          
  ##A2 <- DoAllTheSaving(GFMAA, OneVV, KPAm, PrMeVec, NCount) 
  ##A2b <- DoAllTheSavingL(GFMAA, OneVV, KPAm, PrMeVec, NCount)  
 ##########################
 ### Second Paper ##4
  PrMeVec <- c(4,7,9,8,10); 
  ##GFMAA <- TRAA2;
  ##OneVV <- OneVV;
  NCount <- 200; KPAm <- 4;
   TMM = FALSE; rndit = 2;
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
  ##GFMAA <-   GetForMeAAA2(  OneVV, NumReps = 200, PrintFlag = 0, BRealInput = c(4,3,-2.5,1),  LARSSeekFlag = 2.5 )                                          
  ##A3 <- DoAllTheSaving(GFMAA, OneVV, KPAm, PrMeVec, NCount)
  ##A3b <- DoAllTheSavingL(GFMAA, OneVV, KPAm, PrMeVec, NCount)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (OneVV, NumReps, PrintFlag = 0, BRealInput = c(1, -1, 
    1, -1, 1, -1), LARSSeekFlag = 0, musize = 10) 
{
    ALT2 <- 10 * 11/2
    ALT3 <- 10
    if (PrintFlag < 0) {
        PrintOutFlags <- 0
    }
    else if (PrintFlag == 0) {
        PrintOutFlags <- 1
    }
    else {
        PrintOutFlags <- PrintFlag
    }
    StuffMeMatrix <- matrix(0, length(OneVV[, 1]) * NumReps, 
        7 + ALT2 * 3 + ALT3 * 3)
    SavedOutPutDirectory <- LoadSavedOutPutDirectory()
    FilesInDIR <- unlist(list.files(SavedOutPutDirectory))
    countt <- 0
    for (ii in 1:length(OneVV[, 1])) {
        ppiuse <- OneVV[ii, 1]
        sigmaFactCov <- OneVV[ii, 2]
        ppCov <- OneVV[ii, 3]
        NNOn <- OneVV[ii, 4]
        NROn <- OneVV[ii, 5]
        NoiseOn <- OneVV[ii, 6]
        GMF <- GetForMeFunction2(ppiuse = ppiuse, sigmaFactCov = sigmaFactCov, 
            ppCov = ppCov, NNOn = NNOn, NROn = NROn, NoiseOn = NoiseOn, 
            NumReps = NumReps, FilesInDIR = FilesInDIR, PrintOutFlags = PrintOutFlags, 
            BRealInput = BRealInput, LARSSeekFlag = LARSSeekFlag)
        if (dim(as.matrix(GMF))[2] != length(StuffMeMatrix[1, 
            ])) {
            GMF = as.matrix(GMF)
            print(paste("GetForMeAAA2: Stuff Error at OneVV[", 
                ii, "], we have : ", sep = ""))
            print(paste("Def = 10", "  ALT2 = ", ALT2, "  ALT3 = ", 
                ALT3, sep = ""))
            print(paste(" StuffMeMatrixDims = c( ", length(StuffMeMatrix[, 
                1]), ", ", length(StuffMeMatrix[1, ]), ")", sep = ""))
            print(paste(" GMF = c( ", length(GMF[, 1]), ", ", 
                length(GMF[1, ]), ")", sep = ""))
            print(paste("  NumReps = ", NumReps, sep = ""))
            print(paste("  countt = ", countt, " and countt+NumReps=", 
                countt + NumReps, sep = ""))
            retrn(StuffMeMatrix)
        }
        StuffMeMatrix[(countt + 1):(countt + NumReps), ] <- as.matrix(GMF)
        countt <- countt + NumReps
        if (PrintFlag > 0 && countt\%\%PrintFlag == 0) {
            plot(GMF[, 8] ~ GMF[, 9], main = paste("ppiuse = ", 
                ppiuse, ", sigFC = ", sigmaFactCov, ", ppC = ", 
                ppCov, "\n NNOn = ", NNOn, ", NROn = ", NROn, 
                ", NoiseOn = ", NoiseOn, "countt = ", countt, 
                "\n", "LARSSeekFlag = ", LARSSeekFlag, sep = ""))
        }
    }
    return(StuffMeMatrix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }
\keyword{ datagen }% __ONLY ONE__ keyword per line
