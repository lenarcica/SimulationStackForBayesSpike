\name{MySaveFileName}
\alias{MySaveFileName}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get a File Name to save to }
\description{
  Returns a File Name
}
\usage{
MySaveFileName(OneVV, KPAm, NCount, PrMeVec, LL = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OneVV}{ Matrix of Parameter Input values, there are 5 inputs per sim. }
  \item{KPAm}{ Length of Active True Dataset}
  \item{NCount}{ Number of Simulations }
  \item{PrMeVec}{ What Estimators to Print }
  \item{LL}{ Whether it is a long table or not }
}
\details{
  A File Name for DoAllTheSaving
}
\value{
   A File Name
}
\references{ ~put references to the literature/web site here ~ }
\author{ Alan Lenarcic }
\note{ 
}
\seealso{ DoAllTheSaving, DoAllTheSavingL }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (OneVV, KPAm, NCount, PrMeVec, LL = FALSE) 
{
    STD <- LoadSavedTableDirectory()
    if (LL == TRUE) {
        My = "L"
    }
    else {
        My = ""
    }
    name <- paste(STD, "/", "OutputTable", My, "KP", KPAm, "CNT", 
        NCount, "TB", paste(PrMeVec, collapse = ""), "mNN", tSeq(min(OneVV[, 
            4])), "MNN", tSeq(max(OneVV[, 4])), "mKP", tSeq(min(OneVV[, 
            5])), "MKP", tSeq(max(OneVV[, 5])), "msig", tSeq(min(OneVV[, 
            6])), "Msig", tSeq(max(OneVV[, 6])), ".tex", sep = "")
    return(name)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ documentation }% __ONLY ONE__ keyword per line
