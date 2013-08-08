\name{makecov}
\alias{makecov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make covariance matrix.
}
\description{
Make covariance matrix.
}
\usage{
makecov(theta, correlation, design, subject)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{correlation}{
%%     ~~Describe \code{correlation} here~~
}
  \item{design}{
%%     ~~Describe \code{design} here~~
}
  \item{subject}{
%%     ~~Describe \code{subject} here~~
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
