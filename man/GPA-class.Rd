\name{GPA-class}
\alias{GPA-class}
\alias{show,GPA-method}
\alias{print,GPA-method}
\alias{fdr,GPA-method}
\alias{cov,GPA-method}
\alias{estimates,GPA-method}
\alias{fdr}
\alias{cov}
\alias{estimates}
\Rdversion{1.1}
\docType{class}
\title{Class "GPA" }
\description{
This class represents GPA model fit.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GPA", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"list"}, 
    	representing the fitted GPA model. }
    \item{\code{setting}:}{Object of class \code{"list"}, 
    	representing the setting for GPA model fitting. }
    \item{\code{gwasPval}:}{Object of class \code{"matrix"}, 
    		representing the p-value matrix from GWAS data. }
    \item{\code{annMat}:}{Object of class \code{"matrix"}, 
    		representing the annotation matrix. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "GPA")}: provide brief summary of the object. }
    \item{print}{\code{signature(x = "GPA")}: 
		provide the matrix of posterior probability that a SNP belongs to each combination of association status. }
    \item{fdr}{\code{signature(object = "GPA")}: 
		provide the matrix of local FDR that a SNP is not associated with each phenotype,
		where the order of column is same as those in input GWAS data. }
    \item{cov}{\code{signature(object = "GPA", silent=FALSE, vDigitEst=1000, vDigitSE=1000 )}: 
        provide the covariance matrix for parameter estimates of GPA model.
		If \code{silent=TRUE}, it suppresses the summary output.
		\code{vDigitEst} and \code{vDigitSE} control number of digits for reporting parameter estimates and standard errors. For example, setting it to 1000 means printing out values up to three digits below zero.}
    \item{estimates}{\code{signature(object = "GPA")}: 
        extract estimates from GPA model fit. }
     }
}
\references{
Chung D, Yang C, Li C, Gelernter J, and Zhao H (2014), 
"GPA: A statistical approach to prioritizing GWAS results 
by integrating pleiotropy information and annotation data."
}
\author{ Dongjun Chung }
\seealso{
\code{\link{GPA}}, \code{\link{pTest}}, \code{\link{aTest}}.
}
\examples{
showClass("GPA")

\dontrun{		
fit.GPA.wAnn <- GPA( pmat, ann )
fit.GPA.wAnn
pp.GPA.wAnn <- print( fit.GPA.wAnn )
fdr.GPA.wAnn <- fdr( fit.GPA.wAnn )
cov.GPA.wAnn <- cov( fit.GPA.wAnn )
est.GPA.wAnn <- estimates( fit.GPA.wAnn )
}
}
\keyword{classes}
