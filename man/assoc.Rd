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
\S4method{assoc}{GPA}( object, FDR=0.05, fdrControl="global" )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{GPA model fit. }
    \item{FDR}{FDR level. }
	\item{fdrControl}{Method to control FDR. Possible values are "global" (global FDR control) and
        "local" (local FDR control). Default is "global". } 
    \item{...}{ Other parameters to be passed through to generic \code{assoc}.}
}
\details{
\code{assoc} uses the direct posterior probability approach of Newton et al. (2004)
to control global FDR in association mapping.
}
\value{
Returns a binary matrix indicating association of each SNP,
where its rows and columns match those of input p-value matrix for function \code{GPA}.
}
\references{
Chung D, Yang C, Li C, Gelernter J, and Zhao H (2014), 
"GPA: A statistical approach to prioritizing GWAS results 
by integrating pleiotropy information and annotation data."

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
