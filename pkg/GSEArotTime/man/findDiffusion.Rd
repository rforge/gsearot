\name{findDiffusion}
\alias{findDiffusion}
\title{
Find the diffusion matrix based on a graph.
}
\description{
Find the diffusion matrix based on a graph.
}
\usage{
findDiffusion(graph, alpha = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{
%%     ~~Describe \code{graph} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
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
function (graph, alpha = 1) 
{
    require(igraph)
    graph <- as.undirected(graph)
    A <- get.adjacency(graph)
    degs <- degree(graph)
    L <- A - diag(degs)
    S <- eigen(L)
    K <- S$vectors \%*\% diag(exp(alpha * S$values)) \%*\% t(S$vectors)
    idx <- which(K < 0)
    K[idx] <- 0
    K
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
