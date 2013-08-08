\name{gsea.rotation}
\alias{gsea.rotation}
\title{
GSEA rotation for longitudinal data.
}
\description{
GSEA rotation for longitudinal data. Assumes common covariance matrix for all genes.
}
\usage{
gsea.rotation(S, y, X, contrast, covmat = NULL, nrot = 10000, ES.p = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
Matrix indicating gene set membership with genes as rows and gene sets 
as columns. 1 indicates that the gene is a member, 0 indicates not a member.
}
  \item{y}{
Matrix of gene expression data. Genes are columns and samples are rows.
}
  \item{X}{
Design matrix. Each row represent a sample.
}
  \item{contrast}{
Matrix of contrasts. Each column represent one contrast.
}
  \item{covmat}{
 Covariance matrix. 
}
  \item{nrot}{
number of rotations (default: 10000).
}
  \item{ES.p}{
Weighting parameter for enrichment score (default=1).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{Ngenes}{Size of gene sets.}
\item{ES}{Enrichment scores.}
\item{NES}{Normalised enrichment scores.}
\item{p.value}{p-value per gene set.}
\item{q.value}{q-value per gene set.}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Guro Dorum <guro.dorum@umb.no>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gsea.rotation2}} for GSEA rotation for longitudinal data assuming common covariance 
matrix only within gene sets.
}
\examples{
#2 gene sets of 3 genes each. 3 time points. 1 replicate at each time point. Total number of samples per gene is n
ngs <- 10
ng <- 10
nt <- 3
nrep <- 10
n <- nt*nrep

#Indicating which gene set each gene belongs to
iset <- rep(1:ngs, each=ng)
S <- matrix(0,ng*ngs,ngs)
for(i in 1:ngs) {
  S[iset==i,i] <- 1
}

Timer <- rep(1:3,each=nrep) 
Time <- ordered(Timer)
#Design for random factors
Z <- cbind(Timer)
colnames(Z) <- c("Time")
#Design for fixed factors
X <- model.matrix(~Time)
X <- X[,-ncol(X)]
colnames(X) <- c("Gene","Linear")

#Set values for variance components
sigmat <- 2 #Time variance
sigmae <- 3 #Random error variance
phi <- 0.9 #Time correlation parameter
#Correlation between time points for same gene
Ts <- matrix(0,n,n)
for(i in 1:(n-1)) {
  for(j in (i+1):n) {
      Ts[i,j] <- abs(Z[j,1] - Z[i,1]) 
  }
}
Ts <- Ts + t(Ts)
R <- exp(-phi*Ts)
#Covariance matrix
subject <- rep(1,each=n)
V <- makecov(theta=c(sigmat,sigmae),correlation=R,design=Z,subject=subject)
Vroot <- chol(V)

#Simulate data with covariance structure
Y1 <- matrix(rnorm(n*ng*ngs),ng*ngs,n)
Y <- t(crossprod(t(Y1),Vroot)) 
#Add fixed gene effects and linear time effecs first gene set
betaG <-  2; betaLT <- 1
Yeff <- Y
Yeff[,1:ng] <- Yeff[,1:ng] + matrix(X[,1]*betaG,ncol=ng,nrow=n) + matrix(X[,2]*betaLT,ncol=ng,nrow=n)

#Estimate covariance matrix
#Constraints for optimisation and start values
A <- rbind(c(1,0), c(0,1))
b <- c(1e-6, 1e-6)
alphaStart <- c(sigmae/sigmat,phi)
reml_res <- reml( alphaStart=alphaStart, X=X, Y=Y, Z=Z, subject=subject, A=A, b=b)
sigmat.hat <- reml_res$sigma
sigmae.hat <- reml_res$alpha[1]*sigmat.hat
phi.hat <- reml_res$alpha[2]
R.hat <- exp(-phi.hat*Ts)
#Estimated covariance matrix
V.hat <- makecov(theta=c(sigmat.hat,sigmae.hat),correlation=R.hat,design=Z,subject=subject)

#Contrasts for testing gene effect and linear time effect
contrast <- matrix(0,ncol(X),2)
colnames(contrast) <- c("Gene","Linear")
contrast[colnames(X)=="Gene",1] <- 1
contrast[colnames(X)=="Linear",2] <- 1
      
#GSEA rotation
gsea.rotation(S=S, y=t(Yeff), X=X, contrast=contrast, covmat=V.hat, nrot=1000, ES.p=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
