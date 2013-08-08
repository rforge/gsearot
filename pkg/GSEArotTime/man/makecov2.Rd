\name{makecov2}
\alias{makecov2}
\title{
Make covariance matrix.
}
\description{
Make covariance matrix.
}
\usage{
makecov2(theta, Rt, Rg, design)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{Rt}{
%%     ~~Describe \code{Rt} here~~
}
  \item{Rg}{
%%     ~~Describe \code{Rg} here~~
}
  \item{design}{
%%     ~~Describe \code{design} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
