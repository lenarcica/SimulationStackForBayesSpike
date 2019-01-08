\name{SubRows}
\alias{SubRows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Rows Subfunction }
\description{
  Sub Function that labels row headings for DoAllTheSaving Latex Tables
}
\usage{
SubRows(OneVV, KPAm, rndit = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OneVV}{ OneVV input matrix of vectors}
  \item{KPAm}{ # of active factors for sim }
  \item{rndit}{ # of digits to round to with rd0 }
}
\details{
  DoAllTheSaving
}
\value{
 
}
\references{ Lenarcic Thesis 2009 }
\author{ Alan Lenarcic }
\note{ 
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (OneVV, KPAm, rndit = 2) 
{
    rt <- paste(" \\begin{array}{c}  P_{\\mbox{\\tiny{xcor}}} = ", 
        ".", unlist(strsplit(as.character(round(OneVV[, 2], rndit)), 
            "\\."))[2], " \\mbox{ , } \\xi = ", ".", unlist(strsplit(as.character(round(OneVV[, 
            3], rndit)), "\\."))[2], " \\\\", " \\kappa_{\\mbox{\\tiny{$\\mathcal{A}$}}} = ", 
        KPAm, " \\mbox{ , } \\sigma = ", OneVV[, 6], "\\\\", 
        " n = ", OneVV[, 4], "\\mbox{ , } ", " \\kappa = ", OneVV[, 
            5], "\\end{array} ", sep = "")
    return(rt)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ print }% __ONLY ONE__ keyword per line
