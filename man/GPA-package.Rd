\name{GPA-package}
\alias{GPA-package}
\docType{package}
\title{
GPA (Genetic analysis incorporating Pleiotropy and Annotation)
}
\description{
This package provides functions for fitting GPA,
a statistical approach to prioritizing GWAS results 
by integrating pleiotropy information and annotation data.
}
\details{
\tabular{ll}{
Package: \tab GPA\cr
Type: \tab Package\cr
Version: \tab 0.9.3\cr
Date: \tab 2014-10-31\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
This package contains a main class, \code{GPA}, which represents GPA model fit.
This package contains four main methods, 
\code{GPA}, \code{assoc}, \code{pTest}, and \code{aTest}.
\code{GPA} method fits the GPA model
and \code{assoc} method implements association mapping.
\code{pTest} and \code{aTest} methods implement hypothesis testing
for pleiotropy and annotation enrichment, respectively.
}
\author{
Dongjun Chung

Maintainer: Dongjun Chung <chungd@musc.edu>
}
\references{
Chung D*, Yang C*, Li C, Gelernter J, and Zhao H (2014), 
"GPA: A statistical approach to prioritizing GWAS results 
by integrating pleiotropy information and annotation data."
To appear in PLoS Genetics. (*: Joint first authors)
}
\keyword{ package }
\seealso{
\code{\link{GPA}}, \code{\link{assoc}}, \code{\link{pTest}}, \code{\link{aTest}}, 
\code{\linkS4class{GPA}}.
}
\examples{
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
}