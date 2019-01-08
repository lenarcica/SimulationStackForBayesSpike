\name{DoAllTheSaving}
\alias{DoAllTheSaving}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Save a Latex Table }
\description{
  Saves a Latex Table based upon SimForMe2 Output
}
\usage{
DoAllTheSaving(GFMAA, OneVV, KPAm, PrMeVec, NCount, rndit=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GFMAA}{ Output from SimForMe2 }
  \item{OneVV}{ Matrix of Simulation Parameters used }
  \item{KPAm}{ Active Covariate set used }
  \item{PrMeVec}{ Which estimates, 1-10 to use }
  \item{NCount}{ Count of simulations per parameter set (200 for typical) }
  \item{rndit}{Number of Decimal Places to Round}
}
\details{
  Saves a Latex Table of Simulation output
}
\value{
   None
}
\references{ None }
\author{ Alan Lenarcic }
\note{  Creates a standalone Latex Table
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
  PrMeVec <- c(4,7,9,8,10); 
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ print }% __ONLY ONE__ keyword per line
