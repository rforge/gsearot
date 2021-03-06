\name{reml2}
\alias{reml2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Restricted maximum likelihood to estimate the covariance matrix of longitudinal expression data with gene dependencies within gene sets with correlation between samples.
}
\description{
Restricted maximum likelihood for longitudinal expression data where genes are assumed correlated within gene set and uncorrelated between different gene sets. All genes within the same gene set are assumed to have identical covariance matrix. Also requires that all gene sets are of same size and hence have identical covariance matrix.
}
\usage{
reml2(alphaStart, X, Y, Z, Ts, Gs, Bs = NULL, set, A, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alphaStart}{
A numeric vector of start values for the parameters to be estimated. The first parameter must be gene variance, the second parameter random error variance, the third must be a parameter controlling correlation between time points and the fourth must be a parameter controlling correlation between genes. The remaining parameters can be variances for other random design factors.
}
  \item{X}{
Design matrix with fixed design variables for a gene set (assuming same design for all gene sets)
}
  \item{Y}{
Vector of observed values for all gene sets
}
  \item{Z}{
Design matrix for random design factors for a gene set (assuming same design for all gene sets)
}
  \item{Ts}{
Structure of time dependencies between samples in a gene set (assuming same design for all gene sets)
}
  \item{Gs}{
Structure of gene dependencies between samples in a gene set (assuming same design for all gene sets)
}
  \item{Bs}{
(Optional) structure of batch dependencies between samples in a gene set (assuming same design for all gene sets)
}
  \item{set}{
Vector of indices indicating gene set for each sample
}
  \item{A}{
Restrictions for parameters. See \code{link\constrOptim}}
}
  \item{b}{
Restrictions for parameters. See \code{link\constrOptim}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{alpha.hat}{Vector of estimates for the parameters in alphaStart}
\item{beta.hat}{Vector of estimated parameters for fixed design factors}
\item{sigma}{Estimated common variance}
\item{maklix}{The maximum likelihood value}
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
\code{\link{reml}} for restricted maximum likelihood without gene correlations and \code{\link{reml_genedep}} that does not assume identical gene set size.
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
