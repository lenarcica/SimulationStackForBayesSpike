\name{EMRIDGESolveP}
\alias{EMRIDGESolveP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Old EMRidge Solve P Version }
\description{
  Old Inactive Version of EMRidge that solves for pi_A
}
\usage{
EMRIDGESolveP(yys, xxs, n1 = 128, n2 = 40, EMiters = 10, pstart = -999, p1p = 1, p2p = 2, sigmaalpha = 0.5, sigmatau = 1, tauDsqStart = 1, tauAsqStart = 4, tauDMult = -999, tauAMult = -999, tauAMultStop = 0, BStartBetas = -1, SigmaSqNoiseSt = 0.5, ZeroCriterion = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yys}{ ~~Describe \code{yys} here~~ }
  \item{xxs}{ ~~Describe \code{xxs} here~~ }
  \item{n1}{ ~~Describe \code{n1} here~~ }
  \item{n2}{ ~~Describe \code{n2} here~~ }
  \item{EMiters}{ ~~Describe \code{EMiters} here~~ }
  \item{pstart}{ ~~Describe \code{pstart} here~~ }
  \item{p1p}{ ~~Describe \code{p1p} here~~ }
  \item{p2p}{ ~~Describe \code{p2p} here~~ }
  \item{sigmaalpha}{ ~~Describe \code{sigmaalpha} here~~ }
  \item{sigmatau}{ ~~Describe \code{sigmatau} here~~ }
  \item{tauDsqStart}{ ~~Describe \code{tauDsqStart} here~~ }
  \item{tauAsqStart}{ ~~Describe \code{tauAsqStart} here~~ }
  \item{tauDMult}{ ~~Describe \code{tauDMult} here~~ }
  \item{tauAMult}{ ~~Describe \code{tauAMult} here~~ }
  \item{tauAMultStop}{ ~~Describe \code{tauAMultStop} here~~ }
  \item{BStartBetas}{ ~~Describe \code{BStartBetas} here~~ }
  \item{SigmaSqNoiseSt}{ ~~Describe \code{SigmaSqNoiseSt} here~~ }
  \item{ZeroCriterion}{ ~~Describe \code{ZeroCriterion} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{  Alan Lenarcic }
\note{ 
}
\seealso{ EMRIDGE, EMPI2Lasso }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (yys, xxs, n1 = 128, n2 = 40, EMiters = 10, pstart = -999, 
    p1p = 1, p2p = 2, sigmaalpha = 0.5, sigmatau = 1, tauDsqStart = 1, 
    tauAsqStart = 4, tauDMult = -999, tauAMult = -999, tauAMultStop = 0, 
    BStartBetas = -1, SigmaSqNoiseSt = 0.5, ZeroCriterion = 4) 
{
    if (!is.loaded("EMRIDGEShell2")) {
        LoadTwoCC()
    }
    LoadTwoCC()
    p1p = abs(p1p)
    p2p = abs(p2p)
    if (pstart < 0 || pstart > 1) {
        pstart = p1p/(p1p + p2p)
    }
    puse = pstart
    Returner = list(ReturnBetasPI = matrix(0, EMiters, length(XX[1, 
        ])), BBHatsKeepPI = matrix(0, EMiters, length(XX[1, ])), 
        sigmaNoiseFit = vector("numeric", EMiters), pfit = vector("numeric", 
            EMiters), Elpfit = vector("numeric", EMiters), El1mpfit = vector("numeric", 
            EMiters), sigmaInvNoiseFit = vector("numeric", EMiters), 
        ReturnBetas = ReturnBetas)
    SigmaNoiseForInput = SigmaSqNoiseSt
    E1OverSigma = 1/SigmaNoiseForInput
    Elp = log(puse)
    El1mp = log(1 - puse)
    for (ii in 1:EMiters) {
        MyResult <- EMRIDGE(yys = YY, xxs = XX, TotalRuns = 128, 
            NumEMConv = 40, ppiuse = puse, tauDsqStart = tauDsqStart, 
            tauAsqStart = tauAsqStart, tauDMult = tauDMult, tauAMult = tauAMult, 
            tauAMultStop = tauAMultStop, BStartBetas = -1, SigmaSqNoise = SigmaNoiseForInput, 
            EFlag = TRUE, logp = Elp, log1mp = El1mp, SigmaNoiseInv = E1OverSigma)
        NewNoise <- sum((YY - XX \%*\% MyResult$BEndBetas)^2)
        E1OverSigma <- (sigmaalpha + length(YY) - 1)/(sigmatau + 
            NewNoise)
        SigmaNoiseForInput <- 1/E1OverSigma
        ETrueNoise <- (sigmatau + NewNoise)/(sigmaalpha + length(YY) - 
            3)
        NumNo = length(MyResult$BEndBetas[round(MyResult$BEndBetas, 
            ZeroCriterion) == 0])
        pp <- 1:99999/1e+05
        Elp <- sum(log(pp) * dbeta(pp, p1p + (length(XX[1, ]) - 
            NumNo), p2p + NumNo) * (pp[2] - pp[1]))
        El1mp <- sum(log(1 - pp) * dbeta(pp, p1p + (length(XX[1, 
            ]) - NumNo), p2p + NumNo) * (pp[2] - pp[1]))
        puse <- sum(pp * dbeta(pp, p1p + (length(XX[1, ]) - NumNo), 
            p2p + NumNo) * (pp[2] - pp[1]))
        Returner$ReturnBetasPI[ii, ] <- MyResult$BEndBetas
        Returner$sigmaNoiseFit[ii] <- ETrueNoise
        Returner$pfit[ii] <- puse
        Returner$Elpfit[ii] <- Elp
        Returner$El1mpfit[ii] <- El1mp
        Returner$sigmaInvNoiseFit <- E1OverSigma
        (Returner$BBHatsKeepPI[ii, ] <- xxx) = 1:length(MyResult$ReturnBetasAll[, 
            1, n2])
        yylim = c(min(MyResult$ReturnBetasAll[, , n2]), max(MyResult$ReturnBetasAll[, 
            , n2]))
        kLen = length(XX[1, ])
        cols <- rainbow(kLen)
        plot(x = xxx, y = MyResult$ReturnBetasAll[1, , n2], main = paste("On EMiter = ", 
            ii, sep = ""), type = "n", ylim = yylim)
        for (jjk in 1:kLen) {
            lines(x = xxx, y = MyResult$ReturnBetasAll[jjk, , 
                n2], col = cols[jjk], lty = jjk, lwd = 2)
        }
    }
    Returner$ReturnBetas <- MyResult
    RetCEV$Type = "EMPIRIDGE"
    return(Returner)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
