\name{reml}
\alias{reml}
\title{
Restricted maximum likelihood.
}
\description{
Restricted maximum likelihood.
}
\usage{
reml(alphaStart, X, Y, Z, subject, A, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alphaStart}{
%%     ~~Describe \code{alphaStart} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{subject}{
%%     ~~Describe \code{subject} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
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
function (alphaStart, X, Y, Z, subject, A, b) 
{
    G <- ncol(Y)
    N <- nrow(Y)
    p <- ncol(X)
    res <- constrOptim(theta = alphaStart, f = loglik, grad = NULL, 
        ui = A, ci = b, X = X, Z = Z, Y = Y, subject = subject)
    alpha.hat <- res$par
    maxlik <- res$value
    phi <- alpha.hat[length(alpha.hat)]
    theta.hat <- c(1, alpha.hat[-length(alpha.hat)])
    R <- matrix(0, N, N)
    for (i in 1:(N - 1)) {
        for (j in (i + 1):N) {
            if (subject[i] == subject[j]) 
                R[i, j] <- exp(-phi * (Z[j, 1] - Z[i, 1]))
        }
    }
    R <- R + t(R)
    diag(R) <- 1
    V.hat <- makecov(theta.hat, correlation = R, design = Z, 
        subject = subject)
    Vinv <- solve(V.hat)
    XVinv <- t(X) \%*\% Vinv
    XVinvX <- XVinv \%*\% X
    beta.hat <- solve(XVinvX) \%*\% XVinv \%*\% Y
    r <- Y - X \%*\% beta.hat
    RSS <- sum(diag(t(r) \%*\% Vinv \%*\% (r)))
    sigma.hat <- RSS/(G * (N - p))
    return(list(alpha = alpha.hat, beta = beta.hat, sigma = sigma.hat, 
        maxlik = maxlik))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
