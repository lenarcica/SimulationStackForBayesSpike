\name{rd0}
\alias{rd0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  Rounds numbers appropriately for Latex tables, shrinks lesser digits
}
\usage{
rd0(RoundNumber)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RoundNumber}{ ~A number to round. }
}
\details{
  None Needed, Latex Formatting
}
\value{
 A string for Latex input
}
\references{ "TwoLasso for Bayes Adapted Selection", Lenarcic,
  Journal of Computational and Graphical Statistics, Submitted 2009. }
\author{ Alan Lenarcic~ }
\note{ 
No Notice
}
\seealso{ TwoDigits }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
 rd0(2.53);
 rd0(1);
 rd0(-2.43);
 rd0(200.43);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ documentation }
\keyword{ internal }% __ONLY ONE__ keyword per line
