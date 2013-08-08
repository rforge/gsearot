\name{fdr}
\alias{fdr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
False discovery rate computation for GSEA.
}
\description{
Computes a false discovery rate (FDR) q-value for a gene
set given the observed normalised enrichment scores and a normalised null
distribution. Calculations are done separately for positive and negative
normalised enrichment scores.
}
\usage{
fdr(NES, NESobs, NESnull)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NES}{
Normalised enrichment score for the given gene set.
}
  \item{NESobs}{
Vector of normalised enrichment scores for all gene sets to be tested.
}
  \item{NESnull}{
Matrix of normalised null distribution for all gene sets to be tested.
Each column represent a gene set, and each row represent a permutation.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{q}{FDR q-value for the given gene set.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S.
and Mesirov, J. P (2005) Gene set enrichment analysis: A knowledge-based
approach for interpreting genome-wide expression pro¯les, \emph{PNAS}, \bold{102},
15545-15550.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FDR }
\keyword{ significance }% __ONLY ONE__ keyword per line
