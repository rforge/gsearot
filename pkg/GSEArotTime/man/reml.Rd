\name{reml}
\alias{reml}
\title{
Restricted maximum likelihood to estimate the covariance matrix of longitudinal expression data with correlation between samples.
}
\description{
Restricted maximum likelihood for longitudinal expression data where genes are assumed independent. All genes are assumed to have identical covariance matrix.
}
\usage{
reml(alphaStart, X, Y, Z, subject, A, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alphaStart}{
A numeric vector of start values for the parameters to be estimated. The first parameter must be random error variance, the last must be a parameter controlling correlation between time points. The remaining parameters can be variances for other random design factors 
}
\item{X}{
Design matrix with fixed design variables
}
  \item{Y}{
Matrix of observed values. Rows are samples and columns are variables.
}
  \item{Z}{
Design matrix for random design factors
}
  \item{A}{
Restrictions for parameters. See \code{link\constrOptim}}
}
  \item{b}{
Restrictions for parameters. See \code{link\constrOptim}}
}
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
\code{\link{reml2}} and \code{\link{reml_genedep}} for restricted maximum likelihood where gene correlations are assumed within gene set.
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
