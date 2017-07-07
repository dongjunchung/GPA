\name{GPA-class}
\alias{GPA-class}
\alias{show,GPA-method}
\alias{print,GPA-method}
\alias{fdr,GPA-method}
\alias{cov,GPA-method}
\alias{estimates,GPA-method}
\alias{se,GPA-method}
\alias{fdr}
\alias{cov}
\alias{estimates}
\alias{se}
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
    \item{fdr}{\code{signature(object = "GPA", pattern=NULL)}: provide local FDR.
    By default (i.e., \code{pattern=NULL}),
    it returns a matrix of local FDR that a SNP is not associated with each phenotype (i.e., marginal FDR),
		where the order of columns is same as that in input GWAS data.
    If a pattern is specified, a vector of corresponding local FDR is provided.
    See the details about how users can specify the pattern.
    }
    \item{cov}{\code{signature(object = "GPA", silent=FALSE, vDigitEst=1000, vDigitSE=1000 )}: 
        provide the covariance matrix for parameter estimates of GPA model.
		If \code{silent=TRUE}, it suppresses the summary output.
		\code{vDigitEst} and \code{vDigitSE} control number of digits for reporting parameter estimates and standard errors. For example, setting it to 1000 means printing out values up to three digits below zero.}
    \item{estimates}{\code{signature(object = "GPA")}: 
        extract parameter estimates from GPA model fit. }
    \item{se}{\code{signature(object = "GPA")}: 
        extract standard errors for parameter estimates from GPA model fit. }
     }
}
\details{
When users use \code{fdr} method, users can specify the pattern using 1 and * in \code{pattern} argument,
where 1 and * indicate phenotypes of interest and phenotypes that are not of interest, respectively.
For example, when there are three phenotypes,
\code{pattern="111"} means a SNP associated with all of three phenotypes,
while \code{pattern="11*"} means a SNP associated with the first two phenotypes
(i.e., association with the third phenotype is ignored (averaged out)).
}
\references{
Chung D*, Yang C*, Li C, Gelernter J, and Zhao H (2014), 
"GPA: A statistical approach to prioritizing GWAS results 
by integrating pleiotropy information and annotation data,"
PLoS Genetics, 10: e1004787. (* joint first authors)
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
fdr11.GPA.wAnn <- fdr( fit.GPA.wAnn, pattern="11" )
fdr1..GPA.wAnn <- fdr( fit.GPA.wAnn, pattern="1*" )
cov.GPA.wAnn <- cov( fit.GPA.wAnn )
est.GPA.wAnn <- estimates( fit.GPA.wAnn )
se.GPA.wAnn <- se( fit.GPA.wAnn )
}
}
\keyword{classes}
