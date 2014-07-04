\name{makecov}
\alias{makecov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make covariance matrix.
}
\description{
Makes a covariance matrix for a model with time series effect and other random effetcs. For balanced and unbalanced designs.
}
\usage{
makecov(theta, correlation, design, subject)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
A vector of variances for all random factors, where first element must be time series variance and second element must be the random error variance. The remaining elements can be variances for any other random design factors
}
  \item{correlation}{
A matrix of time correlations between samples (usually zero for samples from different subjects)
}
  \item{design}{
A nxp design matrix for the p random design factors. Must contain time effect in the first column. The time column must be numeric and the time intervals must be of actual size
}
  \item{subject}{
 A numeric vector indicating which subject each sample belongs to
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Covariance matrix.
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
\code{\link{makecov2}} for covariance matrix also including gene dependencies.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (theta, correlation, design, subject) 
{
    N <- nrow(design)
    V <- matrix(0, N, N)
    for (i in 1:(N - 1)) {
        for (j in (i + 1):N) {
            if (subject[i] == subject[j]) {
                V[i, j] <- theta[1] * correlation[i, j]
            }
            if (ncol(design) > 1) {
                for (k in 2:ncol(design)) {
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
