\name{significance}
\alias{significance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Significance level computation in GSEA
}
\description{
Ccomputes a p-value for a gene set.
}
\usage{
significance(ES, ESnull, doplot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ES}{
Observed enrichment score.
}
  \item{ESnull}{
Null distribution.
}
  \item{doplot}{
If \code{TRUE}, a histogram over the normalised null distribution is drawn. Default is \code{FALSE}.
}
}
\details{
Function that computes a p-value for a gene set given an enrichment score 
and a null distribution. Calculations are done separately for positive and
negative enrichment scores. The function also returns normalised observed 
enrichment score and normalised null distribution for the gene set. 
In addition, a histogram over the normalised null distribution can be plotted.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{p.value }{p-value for the gene set.}
 \item{NES }{Normalised observed enrichment score.}
 \item{NESnull }{Normalised estimated null distribution.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Guro D{\o}rum
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--    or do  help(data=index)  for the standard data sets.

## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pvalue }
\keyword{ GSEA }% __ONLY ONE__ keyword per line
