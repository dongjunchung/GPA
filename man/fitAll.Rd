\name{fitAll}
\alias{fitAll}
\title{
Fit GPA model for all possible pairs of GWAS datasets
}
\description{
Fit GPA model and the GPA model under H0 for all possible pairs of GWAS datasets.
}
\usage{
fitAll( pmat, 
  maxIter=2000, stopping="relative", epsStopLL=1e-10, 
  parallel=FALSE, nCore=8 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pmat}{p-value matrix from GWAS data,
		where row and column correspond to SNP and phenotype, respectively. }
	\item{maxIter}{Maximum number of EM iteration. Default is 2000.}  
  \item{stopping}{Stopping rule for EM iteration. 
		Possible values are \code{"absolute"} (based on absolute difference in log likelihood), 
		\code{"relative"} (based on relative difference in log likelihood), or  
		\code{"aitken"} (based on Aitken acceleration-based stopping rule). 
		Default is \code{"relative"}.}
	\item{epsStopLL}{Threshold to stop the EM iteration. Default is 1e-100.} 
	\item{parallel}{Utilize multiple CPUs for parallel computing 
        using \code{"parallel"} package?
        Possible values are \code{TRUE} (utilize multiple CPUs) 
        or \code{FALSE} (do not utilize multiple CPUs).
Default is \code{FALSE} (do not utilize multiple CPUs). } 
	\item{nCore}{Number of CPUs when parallel computing is utilized.} 
}
\details{
\code{fitAll} function fits the GPA model and the GPA model under H0 for all possible pairs of GWAS datasets. Its output can be used as an input for the \code{shinyGPA} function.
}
\value{
A list with 6 elements, including 
  \code{pmat} (original GWAS p-value matrix),
  \code{combs} (a matrix of GWAS pair indices),
  \code{combList} (a matrix of GWAS pair indices),
  \code{pTestPval} (a matrix of pleiotropy test p-values), 
  \code{fitGPA} (a list of the GPA fit for each pair), and
  \code{fitH0} (a list of the GPA fit under H0 for each pair).
}
\references{
Kortemeier E, Ramos PS, Hunt KJ, Kim HJ, Hardiman G, and Chung D (2017),
"ShinyGPA: An interactive and dynamic visualization toolkit for genetic studies."
}
\author{ Dongjun Chung, Emma Kortemeier }
\seealso{
\code{\link{GPA}}, \code{\link{pTest}}, and \code{\link{shinyGPA}}.
}
\examples{
\dontrun{
# simulator function

simulator <- function( risk.ind, nsnp=20000, alpha=0.6 ) {
  
  m <- length(risk.ind)
  
  p.sig <- rbeta( m, alpha, 1 )
  pvec <- runif(nsnp)
  pvec[ risk.ind ] <- p.sig
  
  return(pvec)
}

# run simulation

set.seed(12345)
nsnp <- 10000
alpha <- 0.4
pmat <- matrix( NA, nsnp, 5 )

pmat[,1] <- simulator( c(1:2000), nsnp=nsnp, alpha=alpha )
pmat[,2] <- simulator( c(501:2500), nsnp=nsnp, alpha=alpha )
pmat[,3] <- simulator( c(4001:6000), nsnp=nsnp, alpha=alpha )
pmat[,4] <- simulator( c(4501:7500), nsnp=nsnp, alpha=alpha )
pmat[,5] <- simulator( c(8001:10000), nsnp=nsnp, alpha=alpha )

# Fit GPA for all possible pairs of GWAS datasets

out <- fitAll( pmat )

# Run the ShinyGPA app using the ouput from fitAll()

shinyGPA(out)
}
}
\keyword{models}
\keyword{methods}
