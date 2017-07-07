\name{GPA-package}
\alias{GPA-package}
\docType{package}
\title{
GPA (Genetic analysis incorporating Pleiotropy and Annotation)
}
\description{
This package provides functions for fitting GPA,
a statistical approach to prioritizing GWAS results 
by integrating pleiotropy information and annotation data,
along with ShinyGPA, a visualization toolkit to investigate 
the pleiotropic architecture using GWAS results.
}
\details{
\tabular{ll}{
Package: \tab GPA\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2017-06-30\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
This package contains a main class, \code{GPA}, which represents GPA model fit.

This package contains four main methods for the GPA framework (Chung et al., 2014), 
\code{GPA}, \code{assoc}, \code{pTest}, and \code{aTest}.
\code{GPA} method fits the GPA model
and \code{assoc} method implements association mapping.
\code{pTest} and \code{aTest} methods implement hypothesis testing
for pleiotropy and annotation enrichment, respectively.

This package contains two main methods for the ShinyGPA visualization toolkit (Kortemeier et al., 2017), \code{fitAll} and \code{shinyGPA}. \code{fitAll} function generates all the intermediate results needed to run \code{shinyGPA}. \code{shinyGPA} opens the ShinyGPA interface, which takes the results generated from \code{fitAll} as input.
}
\author{
Dongjun Chung, Emma Kortemeier

Maintainer: Dongjun Chung <chungd@musc.edu>
}
\references{
Chung D*, Yang C*, Li C, Gelernter J, and Zhao H (2014), 
"GPA: A statistical approach to prioritizing GWAS results 
by integrating pleiotropy information and annotation data,"
PLoS Genetics, 10: e1004787. (* joint first authors)

Kortemeier E, Ramos PS, Hunt KJ, Kim HJ, Hardiman G, and Chung D (2017),
"ShinyGPA: An interactive and dynamic visualization toolkit for genetic studies."
}
\keyword{ package }
\seealso{
\code{\link{GPA}}, \code{\link{assoc}}, \code{\link{pTest}}, \code{\link{aTest}}, 
\code{\linkS4class{GPA}}, \code{\link{fitAll}}, \code{\link{shinyGPA}}.
}
\examples{
\dontrun{
library(GPA)

# simulation setting

nBin <- 20000
pi1 <- 0.2
common <- 0.5
betaAlpha <- c( 0.6, 0.6 )
annP <- c( 0.2, 0.4, 0.4, 0.4 )
seed <- 12345

# simulation setting

nCommon <- round( pi1 * common * nBin )
nUniq <- round( pi1 * ( 1 - common ) * nBin )
nBg <- nBin - 2 * nUniq - nCommon

# M * K matrix of GWAS p-value

set.seed( seed )

pvec1 <- c( rbeta( nCommon, betaAlpha[1], 1 ), rbeta( nUniq, betaAlpha[1], 1 ), 
	runif( nUniq ), runif( nBg ) )
pvec2 <- c( rbeta( nCommon, betaAlpha[2], 1 ), runif( nUniq ),
	rbeta( nUniq, betaAlpha[2], 1 ), runif( nBg ) )
pmat <- cbind( pvec1, pvec2 )

# M * D matrix of annotation
 	
ann <- c( 
	sample( c(1,0), nCommon, replace=TRUE, prob = c( annP[4], 1 - annP[4] ) ), 
	sample( c(1,0), nUniq, replace=TRUE, prob = c( annP[2], 1 - annP[2] ) ),
	sample( c(1,0), nUniq, replace=TRUE, prob = c( annP[3], 1 - annP[3] ) ),
	sample( c(1,0), nBg, replace=TRUE, prob = c( annP[1], 1 - annP[1] ) ) )
		
# GPA without annotation data

fit.GPA.noAnn <- GPA( pmat, NULL )
cov.GPA.noAnn <- cov( fit.GPA.noAnn )
		
# GPA with annotation data

fit.GPA.wAnn <- GPA( pmat, ann )
cov.GPA.wAnn <- cov( fit.GPA.wAnn )

# GPA under pleiotropy H0

fit.GPA.pleiotropy.H0 <- GPA( pmat, NULL, pleiotropyH0=TRUE )

# association mapping

assoc.GPA.wAnn <- assoc( fit.GPA.wAnn, FDR=0.05, fdrControl="global" )

# hypothesis testing for pleiotropy

test.pleiotropy <- pTest( fit.GPA.noAnn, fit.GPA.pleiotropy.H0 )

# hypothesis testing for annotation enrichment

test.annotation <- aTest( fit.GPA.noAnn, fit.GPA.wAnn )

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