\name{reml_genedep}
\alias{reml_genedep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Restricted maximum likelihood to estimate the covariance matrix of longitudinal expression data with gene dependencies within gene sets and correlation between samples.
}
\description{
Restricted maximum likelihood for longitudinal expression data where genes are assumed correlated within gene set and uncorrelated between different gene sets. All genes within the same gene set are assumed to have identical covariance matrix. The gene sets may differ in size and hence have different covariance matrix.
}
\usage{
reml_genedep(alphaStart, Xlist, Ylist, Zlist, TsList, GsList, BsList, A, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alphaStart}{
A numeric vector of start values for the parameters to be estimated. The first parameter must be gene variance, the second parameter random error variance, the third must be a parameter controlling correlation between time points and the fourth must be a parameter controlling correlation between genes. The remaining parameters can be variances for other random design factors.
}
  \item{Xlist}{
Design matrix with fixed design variables for a gene set (assuming same design for all gene sets)
}
  \item{Ylist}{
Vector of observed values for all gene sets
}
  \item{Zlist}{
Design matrix for random design factors for a gene set (assuming same design for all gene sets)
}
  \item{TsList}{
Structure of time dependencies between samples in a gene set (assuming same design for all gene sets)
}
  \item{GsList}{
Structure of gene dependencies between samples in a gene set (assuming same design for all gene sets)
}
  \item{BsList}{
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
\code{\link{reml}} for restricted maximum likelihood without gene correlations and \code{\link{reml2}} for a function with simpler input format if gene sets are of identical size.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (alphaStart, Xlist, Ylist, Zlist, TsList, GsList, BsList, 
    A, b) 
{
    p <- ncol(Xlist[[1]])
    nset <- length(Xlist)
    res <- constrOptim(theta = alphaStart, f = loglik_genedep, 
        grad = NULL, ui = A, ci = b, X = Xlist, Z = Zlist, Y = Ylist, 
        Ts = TsList, Gs = GsList, Bs = BsList)
    alpha.hat <- res$par
    maxlik <- res$value
    phi.hat <- alpha.hat[3]
    ga.hat <- alpha.hat[4]
    theta.hat <- c(1, alpha.hat[-c(3, 4)])
    RSS <- numeric(nset)
    N <- 0
    for (k in 1:nset) {
        Y <- Ylist[[k]]
        X <- Xlist[[k]]
        Z <- Zlist[[k]]
        Bs <- BsList[[k]]
        Gs <- GsList[[k]]
        Ts <- TsList[[k]]
        n <- length(Y)
        dim(Y) <- c(n, 1)
        N <- N + n
        Rt <- exp(-phi.hat * Ts)
        Rg <- exp(-ga.hat * Gs)
        if (ncol(Z) > 2) {
            V.hat <- theta.hat[4] * Bs + theta.hat[1] * Rt + 
                theta.hat[2] * Rg + diag(n) * theta.hat[3]
        }
        else V.hat <- theta.hat[1] * Rt + theta.hat[2] * Rg + 
            diag(n) * theta.hat[3]
        Vinv <- solve(V.hat)
        XVinv <- t(X) \%*\% Vinv
        XVinvX <- XVinv \%*\% X
        XV <- solve(XVinvX) \%*\% XVinv
        beta.hat <- XV \%*\% Y
        r <- Y - X \%*\% beta.hat
        RSS[k] <- t(r) \%*\% Vinv \%*\% r
    }
    sigma.hat <- sum(RSS)/(N - nset * p)
    return(list(alpha = alpha.hat, beta = beta.hat, sigma = sigma.hat, 
        maxlik = maxlik))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
