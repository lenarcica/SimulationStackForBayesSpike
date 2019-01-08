\name{PenaltyFunction}
\alias{PenaltyFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  A unused Lasso Penalty Function
}
\usage{
PenaltyFunction(yys, Betas, xxs, v1, v2, BetOlds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

}
\details{~~ If necessary, more details than the description above ~~
}
\value{

}
\references{ ~put references to the literature/web site here ~ }
\author{ Alan Lenarcic }
\note{ 
}
\seealso{  }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (yys, Betas, xxs, v1, v2, BetOlds) 
{
    sum((yys - xxs \%*\% Betas)^2) + v1 * sum(abs(Betas)) + v2 * 
        sum(abs(Betas - BetOlds))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ methods }% __ONLY ONE__ keyword per line
