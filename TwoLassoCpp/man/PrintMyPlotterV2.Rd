\name{PrintMyPlotterV2}
\alias{PrintMyPlotterV2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to Plot Simulation Output }
\description{
  Plots Type I, Type II error Plot of Errors
}
\usage{
PrintMyPlotterV2(GFMAA,  OneVV, NumReps = 100, PrMeVec, TMM = FALSE, colList, AlphList, xlims = -999, ylims = -999)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GFMAA}{ Output from GetForMeAAA2 }
  \item{OneVV}{ OneVV set of interesting parameters }
  \item{NumReps}{ Number of repetitions used in simulations }
  \item{PrMeVec}{ Which estimators 2-11 to use }
  \item{TMM}{ Flag of that decides what to do with missing/failed estimates }
  \item{colList}{ List of colors to use of dots }
  \item{AlphList}{ List of letters to use for labels }
  \item{xlims}{ Limits for x axis }
  \item{ylims}{ Limits  for y axis }
}
\details{
  Not really
}
\value{
  Nothing
}
\references{ TwoLasso For Bayes Adapted Selection }
\author{ Alan Lenarcic }
\note{  Helps Plotting
}
\seealso{ PrintMyPlotter, PrintHists2 }
\examples{
  OneVV <- rbind( c(7, .2, .3, 25, 10, .5), 
                       c(7, .2, .3, 25, 20, .5),
                       c(7, .2, .3, 25, 30, .5),
                       c(7, .2, .3, 25, 40, .5),
                       c(7, .2, .3, 15, 50, .5),
                       c(7, .2, .3, 25, 50, .5),
                       c(7, .2, .3, 35, 50, .5),
                       c(7, .2, .3, 45, 50, .5),
                       c(7, .2, .3, 25, 100, .5),
                       c(7, .2, .3, 50, 100, .5),                       
                       c(7, .2, .3, 75, 100, .5),
                       c(7, .2, .3, 100, 100, .5))
##TRAA2 <-   GetForMeAAA2(  OneVV, NumReps = 200, PrintFlag = 0, BRealInput = c(1,-1,1,-1,1,-1) )            
 ##MyV <- GiveMeIDOnOutput(GFMAA = TRAA2, OneVV = OneVV, VNum = 1, PrMeNum= 9);
 colPermaList <- c("red", "pink", "purple", "orange", "blue", "green", "cyan", 
         "magenta", "lightblue",
         "aquamarine")
 AlphPermaList <- c("A", "B", "C", "D", "E", "F", "G", "I", "H", "J");         
      PrMeVec <- c(1,2,3,4,5,6,7,9, 8, 10);
      colList <- colPermaList[PrMeVec];
      AlphList <- AlphPermaList[PrMeVec];
      xlims <- c(0,6); ylims = c(0,6);   
      par(plt = c(.1, .96, .15, .8))
 ##     PrintMyPlotterV2(TRAA2,  OneVV, NumReps = 200, PrMeVec, TMM = FALSE,
 ##    colList, AlphList, xlims = xlims, ylims = -999)   

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
