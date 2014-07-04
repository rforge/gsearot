\name{makecov2}
\alias{makecov2}
\title{
Make covariance matrix.
}
\description{
Make covariance matrix for model with time series correlation, gene correlation and correlations due to other random effetcs. For balanced and unbalanced designs.
}
\usage{
makecov2(theta, Rt, Rg, design)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
A vector of variances for all random factors, where first element must be time series variance, second element must be gene variance, and third element must be the random error variance. The remaining elements can be variances for any other random design factors
}
  \item{Rt}{
A matrix of time correlations between samples (usually zero for samples from different subjects)
}
  \item{Rg}{
A matrix of correlations between genes
}
  \item{design}{
A nxp design matrix for the p random design factors. Must contain time in the first column and gene in the second. The time column must be numeric and the time intervals must be of actual size.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Covariance matrix
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
\code{\link{makecov}} for covariance matrix without gene dependencies.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (theta, Rt, Rg, design) 
{
    N <- nrow(design)
    V <- matrix(0, N, N)
    for (i in 1:(N - 1)) {
        for (j in (i + 1):N) {
            V[i, j] <- theta[1] * Rt[i, j] + theta[2] * Rg[i, 
                j]
            if (ncol(design) > 2) {
                for (k in 3:ncol(design)) {
                  if (design[i, k] == design[j, k]) {
                    V[i, j] <- V[i, j] + theta[k + 1]
                  }
                }
            }
        }
    }
    V <- V + t(V)
    diag(V) <- sum(theta)
    V
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
