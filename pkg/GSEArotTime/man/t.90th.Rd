\name{t.90th}
\alias{t.90th}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes moderated gene wise t-statistics.
}
\description{
Computes gene wise t-statistic with variance stabilised by adding the 90th percentile.
}
\usage{
t.90th(X, y = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Gene expression matrix with samples as rows and genes as columns
}
  \item{y}{
Optional vector of 0/1 indicating sample phenotype
}
}
\details{
Computes a gene wise t-statistic with variance stabilised by adding the 90th percentile 
of all the genes' standard deviations to a gene's standard devitation (Efron et al. 2001).
If y is given, a two-sample test is performed, otherwise a one-sample test is performed.
For a two-sample t-test the approach of Pan et al. (2003) is used.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{tS }{Vector with one smoothed t-value for each gene}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{

Efron, B., Tibshirani, R., Storey, J. D. and Tusher, V. (2001) 
Empirical Bayes analysis of a microarray experiment, 
\emph{J Amer Statist Assoc}, 
\bold{96}, 1151-1160. 

Pan,W., Lin,J. and Le,C.T. (2003)
A mixture model approach to detecting differentially expressed genes with microarray data, 
\emph{Funct Integr Genomics}, 
\bold{3}, 117?124.
}
\author{
Guro Dorum
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
##--    or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Test-statistic }
\keyword{ Smoothed }% __ONLY ONE__ keyword per line
