\name{GetForMeAAA}
\alias{GetForMeAAA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  Ole Method for simulating Data
}
\usage{
GetForMeAAA(OneVV, NumReps, PrintFlag = 0, musize = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OneVV}{ OneVV input matrix }
  \item{NumReps}{ ~~Describe \code{NumReps} here~~ }
  \item{PrintFlag}{ ~~Describe \code{PrintFlag} here~~ }
  \item{musize}{ ~~Describe \code{musize} here~~ }
}
\details{
  Old version for simulating data as specified in OneVV matrix
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{  }
\author{ Alan Lenarcic }
\note{ Not Really
}
\seealso{ GetForMeAAA2, GetForMeAll2 }
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

## The function is currently defined as
function (OneVV, NumReps, PrintFlag = 0, musize = 10) 
{


    StuffMeMatrix <- matrix(0, length(OneVV[, 1]) * NumReps, 
        7 + ALT * 3 + 1)
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
        GMF <- GetForMeFunction(ppiuse, sigmaFactCov, ppCov, 
            NNOn, NROn, NoiseOn, NumReps, FilesInDIR)
        StuffMeMatrix[(countt + 1):(countt + NumReps), ] <- as.matrix(GMF)
        countt <- countt + NumReps
        if (PrintFlag > 0 && countt\%\%PrintFlag == 0) {
            plot(GMF[, 8] ~ GMF[, 9], main = paste("ppiuse = ", 
                ppiuse, ", sigFC = ", sigmaFactCov, ", ppC = ", 
                ppCov, "\n NNOn = ", NNOn, ", NROn = ", NROn, 
                ", NoiseOn = ", NoiseOn, "countt = ", countt, 
                sep = ""))
        }
    }
    return(StuffMeMatrix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ methods }% __ONLY ONE__ keyword per line
