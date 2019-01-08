\name{CoordinateDescentPrototype}
\alias{CoordinateDescentPrototype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Coordinate Descent R code Prototype }
\description{
 A function fitting CoordinateDescent Lasso using native R code.
}
\usage{
CoordinateDescentPrototype(xx = -1, yy = -1, XTX = -1, XTY = -1, NLen = -1, 
                    TotalLoops = 20, MaxEpsilon = .0001, OnBeta = -1, 
                    OnGamma = 1, RecordBetasFlag = FALSE, OnGammas = -999)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{ Covariate set  NLen * klen matrix }
  \item{yy}{ Data set, NLen * 1 vector }
  \item{XTX}{ Covariate Matrix \eqn{\mathbf{X}^{T} \mathbf{X}}{transpose of xx times xx} if smaller, \eqn{klen \times klen}{klen * klen} vector }
  \item{XTY}{ Correlated residuals \eqn{\mathbf{X}^{T}\mathbf{X}}{transpose of xx times yy} is a \eqn{\kappa \times 1}{klen * 1} vector }
  \item{NLen}{ NLen, length of yy if XTY is input}
  \item{TotalLoops}{ Maximum Total Loops allowed by algorithm }
  \item{MaxEpsilon}{ Maximum Total Epsilon per loop to justify convergence }
  \item{OnBeta}{ Input Beta if desired }
  \item{OnGamma}{ A single Gamma value for penalty }
  \item{RecordBetasFlag}{ TRUE if larger record should be returned }
  \item{OnGammas}{ A set of Gammas for each of k covariates }
}
\details{
  Runs Friedman's Coordinate Descent Lasso
}
\value{
    Returns List of Many Values
    \item{ReturnBetas}{Result vector}
    \item{FinalBeta}{Also final result vector}
    \item{XTX}{t(xx) \%*\% xx matrix}
    \item{XTY}{t(xx) \%*\% yy vector}
    \item{OnRecordBetas}{Matrix delineating output from all steps of ccordinate
           descent, useful for checking that TotalRuns/NumCDOConv parameter
           is large enough}    
}
\references{ Friedman 2006 }
\author{ Alan Lenarcic }
\note{ Please input xx and yy or XTX, XTY and NLen
}
\seealso{ CoordinateDescent }
\examples{
library(lars)
data(diabetes)

CDO <- CoordinateDescentPrototype(xx=diabetes$x,
     yy=diabetes$y,
     OnGammas = rep(3, length(diabetes$x[1,]))
      );
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (xx = -1, yy = -1, XTX = -1, XTY = -1, NLen = -1, TotalLoops = 20, 
    MaxEpsilon = 1e-04, OnBeta = -1, OnGamma = 1, RecordBetasFlag = FALSE, 
    OnGammas = -999) 
{
    XTXFlag = -1
    if (is.matrix(xx) == FALSE || length(xx) == 1 || xx == -1) {
        if (is.matrix(XTX) == FALSE || length(XTX) == 1 || XTX == 
            -1) {
            print("CoordinateDescent, Error no proper X given, XTX False")
            return
        }
        if ((is.matrix(XTY) == FALSE && is.vector(XTY) == FALSE) || 
            length(XTY) == 1 || XTY == -1) {
            print("CoordinateDescent, Error no proper X given, XTY False")
            return
        }
        if (length(NLen) > 1 || NLen == -1) {
            print("CoordinateDescent, Error, no proper NLen given")
            return
        }
        XTXFlag = 1
        XXX = XTX
        YYY = XTY
        TNLen = -NLen
        kLen = length(XTX[1, ])
    }
    else if (is.matrix(yy) == FALSE || length(yy) == 1 || yy == 
        -1) {
        print("CoordinateDescent, Error no proper Y given")
        return
    }
    else {
        XTXFlag = 0
        XXX = xx
        YYY = yy
        XTX <- t(xx) \%*\% xx
        XTY <- t(xx) \%*\% yy
        TNLen = NLen
        kLen = length(xx[1, ])
    }
    if (length(OnBeta) != kLen) {
        OnBeta = rep(0, kLen)
    }
    if (length(OnGammas) == 1 || OnGammas == -1) {
        OnGammas = rep(1, kLen)
        OnGammas[1] = -999
        OnGammas = rep(OnGamma, kLen)
    }
    if (RecordBetasFlag == TRUE) {
        OnRecordBetas = matrix(0, TotalLoops, kLen)
    }
    else {
        OnRecordBetas = -999
    }
    if (!is.matrix(XTX) || XTX == -1) {
        XTX <- t(xx) \%*\% xx
        XTY <- t(xx) \%*\% yy
    }
    XTYResid <- XTY
    OnBetas = rep(0, kLen)
    PrevBetas = OnBetas
    for (ii in 1:TotalLoops) {
        PrevBetas = OnBetas
        for (jj in 1:kLen) {
            PRD = XTYResid[jj] + XTX[jj, jj] * OnBetas[jj]
            PropBeta = 0
            if (PRD > OnGammas[jj]) {
                PropBeta = (PRD - OnGammas[jj])/XTX[jj, jj]
            }
            else if (PRD < -OnGammas[jj]) {
                PropBeta = (PRD + OnGammas[jj])/XTX[jj, jj]
            }
            else {
                PropBeta = 0
            }
            XTYResid = XTYResid - XTX[, jj] * (PropBeta - OnBetas[jj])
            OnBetas[jj] = PropBeta
        }
        if (is.matrix(OnRecordBetas)) {
            OnRecordBetas[ii, ] <- OnBetas
        }
        if (sum(abs(OnBetas - PrevBetas)) <= MaxEpsilon) {
            break
        }
    }
    CDO <- list(NLen = NLen, kLen = kLen, xx = xx, yy = yy, XTX = XTX, 
        XTY = XTY, OnBeta = OnBetas, OnGamma = OnGamma, OnRecordBetas = OnRecordBetas, 
        FinalBeta = OnBetas, TotalLoops = 20, MaxEpsilon = 3e-04, 
        OnGammas = OnGammas, ReturnBetas = OnBetas, NumLoops = ii)
    return(CDO)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
