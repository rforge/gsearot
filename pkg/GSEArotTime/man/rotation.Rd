\name{rotation}
\alias{rotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rotations
}
\description{
Rotates a gene expression matrix as a step for generating a null distribution for a given gene set.
}
\usage{
rotation(Y, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Matrix with expression values, where rows are samples and columns are genes.
}
  \item{method}{
The method to be used for rotation, either 1 (all rotated genes have expectation 0)
or 2 (non-zero means allowed)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{X }{Matrix with rotated expression values.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Guro Dorum
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Rotation }
\keyword{ GSEA }% __ONLY ONE__ keyword per line
