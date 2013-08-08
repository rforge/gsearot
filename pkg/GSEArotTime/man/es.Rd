\name{es}
\alias{es}
\title{
Enrichment score
}
\description{
Calculates an enrichment score for a gene set.
}
\usage{
es(L, S, p = 1, doplot = FALSE, index = FALSE)
}
\arguments{
  \item{L}{
Vector with one expression value per gene.
}
  \item{S}{
Vector of same length as L where 1 indicates that the gene in L is present in the gene set and 1 indicates that it is not.
}
  \item{p}{
Weight. Default is 1.
}
  \item{doplot}{
If TRUE, the running sum is plotted. Default is FALSE.
}
  \item{index}{
If TRUE, the index at which the ES occurs in the sorted list is returned. Default is FALSE.
}
}
\details{
See Subramanian et al. for details.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{E}{Enrichment score.}
  \item{ind }{Index of the enrichment score (if \code{index=TRUE}).}
%% ...
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S.
and Mesirov, J. P (2005) Gene set enrichment analysis: A knowledge-based
approach for interpreting genome-wide expression pro?les, \emph{PNAS}, \bold{102},
15545-15550.
}
\author{
Guro D\u00F8 rum
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
\keyword{ Enrichment }
\keyword{ Score }% __ONLY ONE__ keyword per line
