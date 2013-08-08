\name{gsea.rotation}
\alias{gsea.rotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GSEA with rotation test
}
\description{
GSEA with rotation test
}
\usage{
gsea.rotation(X, y = NULL, S, nrot = 500, fun = t.90th)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Gene expression matrix with samples as rows and genes as columns.
}
  \item{y}{
Optional vector of 0/1 indicating sample phenotype.
}
  \item{S}{
Matrix indicating gene set membership with genes as rows and gene sets 
as columns. 1 indicates that the gene is a member, 0 indicates not a member.
}
  \item{nrot}{
The number of rotations to perform. Default is 500.
}
  \item{fun}{
Name of function for computing a gene wise test statistic.
Default is t-values with 90th percentile stabilised variance.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{p.values}{A p-value for each gene set}
  \item{q.values}{A FDR q-value for each gene set}
%% ...
}
\references{

D{\o}rum, G., Snipen, L., Solheim, M. and S{\ae}b{\o} (2009) Rotation Testing in Gene Set Enrichment
Analysis for Small Direct Comparison Experiments. \emph{Statistical Applications in Genetics
and Molecular Biology}, \bold{8}(1), article 34.

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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GSEA }
\keyword{ Rotation }% __ONLY ONE__ keyword per line
