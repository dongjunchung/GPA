\name{assoc}
\alias{assoc}
\alias{assoc,GPA-method}
\title{
Association mapping
}
\description{
Association mapping.
}
\usage{
assoc( object, ... )
\S4method{assoc}{GPA}( object, FDR=0.05, fdrControl="global", pattern=NULL )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{GPA model fit. }
    \item{FDR}{FDR level. }
	\item{fdrControl}{Method to control FDR. Possible values are "global" (global FDR control) and
        "local" (local FDR control). Default is "global". } 
    \item{pattern}{Pattern for association mapping.
      By default (i.e., \code{pattern=NULL}),
      \code{assoc} returns a binary matrix indicating association of SNPs for each phenotypes.
      If a pattern is specified, a corresponding binary vector is provided.
      See the details about how users can specify the pattern.}
    \item{...}{ Other parameters to be passed through to generic \code{assoc}.}
}
\details{
\code{assoc} uses the direct posterior probability approach of Newton et al. (2004)
to control global FDR in association mapping.

Users can specify the pattern using 1 and * in \code{pattern} argument,
where 1 and * indicate phenotypes of interest and phenotypes that are not of interest, respectively.
For example, when there are three phenotypes,
\code{pattern="111"} means a SNP associated with all of three phenotypes,
while \code{pattern="11*"} means a SNP associated with the first two phenotypes
(i.e., association with the third phenotype is ignored (averaged out)).
}
\value{
If \code{pattern=NULL}, returns a binary matrix indicating association of SNPs for each phenotype,
where its rows and columns match those of input p-value matrix for function \code{GPA}.
Otherwise, returns a binary vector indicating association of SNPs for the phenotype combination of interest.
}
\references{
Chung D*, Yang C*, Li C, Gelernter J, and Zhao H (2014), 
"GPA: A statistical approach to prioritizing GWAS results 
by integrating pleiotropy information and annotation data,"
PLoS Genetics, 10: e1004787. (* joint first authors)

Newton MA, Noueiry A, Sarkar D, and Ahlquist P (2004),
"Detecting differential gene expression with a semiparametric hierarchical mixture method,"
\emph{Biostatistics}, Vol. 5, pp. 155-176.
}
\author{ Dongjun Chung }
\seealso{
\code{\link{GPA}}, 
\code{\linkS4class{GPA}}.
}
\examples{
\dontrun{
fit.GPA.wAnn <- GPA( pmat, ann )
cov.GPA.wAnn <- cov( fit.GPA.wAnn )
assoc.GPA.wAnn <- assoc( fit.GPA.wAnn, FDR=0.05, fdrControl="global" )
}
}
\keyword{models}
\keyword{methods}
