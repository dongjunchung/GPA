\name{GPA}
\alias{GPA}
\title{
Fit GPA model
}
\description{
Fit GPA model.
}
\usage{
GPA( gwasPval, annMat=NULL, pleiotropyH0=FALSE, empiricalNull=FALSE,
	maxIter=2000, stopping="relative", epsStopLL=1e-10,
	initBetaAlpha=0.1, initPi=0.1, initQ=0.75,
	lbPi=NA, lbBetaAlpha=0.001, lbQ=0.001, lbPval=1e-30,
	vDigit=1000, verbose=1 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{gwasPval}{p-value matrix from GWAS data,
		where row and column correspond to SNP and phenotype, respectively. }
    \item{annMat}{Binary matrix from annotation data,
		where row and column correspond to SNP and annotation, respectively. }
    \item{pleiotropyH0}{Fit GPA under the null hypothesis of pleiotropy test?
        Possible values are \code{TRUE} (under the null hypothesis of pleiotropy test) 
        or \code{FALSE} (usual assumption for GPA).
        Default is \code{FALSE}. }        
    \item{empiricalNull}{Empirically estimate null distribution for GPA?
        Possible values are \code{TRUE} (empirical null distribution) 
        or \code{FALSE} (theoretical null distribution).
        Default is \code{FALSE}. }        
	\item{maxIter}{Maximum number of EM iteration. Default is 2000.}  
    \item{stopping}{Stopping rule for EM iteration. 
		Possible values are \code{"absolute"} (based on absolute difference in log likelihood), 
		\code{"relative"} (based on relative difference in log likelihood), or  
		\code{"aitken"} (based on Aitken acceleration-based stopping rule). 
		Default is \code{"relative"}.}
	\item{epsStopLL}{Threshold to stop the EM iteration. Default is 1e-100.}  
	\item{initBetaAlpha}{Initial value for alpha estimate. Default is 0.1.}  
	\item{initPi}{Initial value for pi estimate. Default is 0.1.} 
	\item{initQ}{Initial value for q estimate. Default is 0.75.} 
	\item{lbPi}{Lower bound for pi estimate. 
		If \code{lbPi=NA}, lower bound is set to 1 / [number of SNPs].
		Default is \code{NA}.}  
	\item{lbBetaAlpha}{Lower bound for alpha estimate. Default is 0.001.}  
	\item{lbQ}{Lower bound for q estimate. Default is 0.001.}  
	\item{lbPval}{Lower bound for GWAS p-value. 
		Any GWAS p-values smaller than \code{lbPval} are set to \code{lbPval}.
		Default is 1e-30.}  
	\item{vDigit}{Number of digits for reporting parameter estimates.
		For example, setting it to 1000 means printing out values up to three digits below zero.} 
	\item{verbose}{Amount of progress report during the fitting procedure. 
		Possible values are 0 (minimal output), 1, 2, or 3 (maximal output).
		Default is 1.} 
}
\details{
\code{GPA} fits the GPA model. It requires to provide GWAS p-value to \code{gwasPval},
while users can also provide annotation data to \code{annMat}.
It is assumed that number of rows of matrix provided to \code{gwasPval} 
equals to that provided to \code{annMat}.

\code{pTest} implements the hypothesis testing for pleiotropy.
It requires two GPA model fits, one of interest and one under the null hypothesis,
and they can be obtained by setting \code{pleiotropyH0=FALSE} and \code{pleiotropyH0=TRUE},
respectively.

\code{aTest} implements the hypothesis testing for annotation enrichment.
It requires two GPA model fits, 
one fitted with using annotation data and one fitted without using annotation data,
and they can be obtained by providing annotation data to \code{annMat} and not, respectively.
}
\value{
Construct \code{GPA} class object.
}
\references{
Chung D*, Yang C*, Li C, Gelernter J, and Zhao H (2014), 
"GPA: A statistical approach to prioritizing GWAS results 
by integrating pleiotropy information and annotation data,"
PLoS Genetics, 10: e1004787. (* joint first authors)
}
\author{ Dongjun Chung }
\seealso{
\code{\link{assoc}}, \code{\link{pTest}}, \code{\link{aTest}}, 
\code{\linkS4class{GPA}}.
}
\examples{
\dontrun{
# GPA without annotation data

fit.GPA.noAnn <- GPA( pmat, NULL )
cov.GPA.noAnn <- cov( fit.GPA.noAnn )
		
# GPA with annotation data

fit.GPA.wAnn <- GPA( pmat, ann )
cov.GPA.wAnn <- cov( fit.GPA.wAnn )

# GPA under the null hypothesis of pleiotropy test

fit.GPA.pleiotropy.H0 <- GPA( pmat, NULL, pleiotropyH0=TRUE )
}
}
\keyword{models}
\keyword{methods}
