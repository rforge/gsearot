\name{reml2}
\alias{reml2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Restricted maximum likelihood.
}
\description{
Restricted maximum likelihood.
}
\usage{
reml2(alphaStart, X, Y, Z, Ts, Gs, Bs = NULL, set, A, b)
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
  \item{Ts}{
%%     ~~Describe \code{Ts} here~~
}
  \item{Gs}{
%%     ~~Describe \code{Gs} here~~
}
  \item{Bs}{
%%     ~~Describe \code{Bs} here~~
}
  \item{set}{
%%     ~~Describe \code{set} here~~
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
function (alphaStart, X, Y, Z, Ts, Gs, Bs = NULL, set, A, b) 
{
    N <- length(Y)
    p <- ncol(X)
    nset <- length(unique(set))
    n <- N/nset
    if (is.null(Bs)) 
        Bs <- matrix(0, ncol(Ts), ncol(Gs))
    res <- constrOptim(theta = alphaStart, f = loglik2, grad = NULL, 
        ui = A, ci = b, X = X, Z = Z, Y = Y, Ts = Ts, Gs = Gs, 
        Bs = Bs, set = set)
    alpha.hat <- res$par
    maxlik <- res$value
    phi <- alpha.hat[3]
    ga <- alpha.hat[4]
    theta.hat <- c(1, alpha.hat[-c(3, 4)])
    Rt <- exp(-phi * Ts)
    Rg <- exp(-ga * Gs)
    if (ncol(Z) > 2) {
        V.hat <- theta.hat[4] * Bs + theta.hat[1] * Rt + theta.hat[2] * 
            Rg + diag(n) * theta.hat[3]
    }
    else V.hat <- makecov2(theta.hat, Rt = Rt, Rg = Rg, design = Z)
    Vinv <- solve(V.hat)
    XVinv <- t(X) \%*\% Vinv
    XVinvX <- XVinv \%*\% X
    XV <- solve(XVinvX) \%*\% XVinv
    RSS <- numeric(nset)
    for (k in 1:nset) {
        idx <- which(set == k)
        Yk <- Y[idx]
        beta.hat <- XV \%*\% Yk
        r <- Yk - X \%*\% beta.hat
        RSS[k] <- t(r) \%*\% Vinv \%*\% r
    }
    sigma.hat <- sum(RSS)/(nset * (n - p))
    return(list(alpha = alpha.hat, beta = beta.hat, sigma = sigma.hat, 
        maxlik = maxlik))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
