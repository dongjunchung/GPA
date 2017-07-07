\name{shinyGPA}
\alias{shinyGPA}
\title{
Run ShinyGPA app
}
\description{
Run ShinyGPA app.
}
\usage{
shinyGPA( out=NULL )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{out}{output of \code{fitAll} function. }
}
\details{
\code{shinyGPA} runs the ShinyGPA app. It takes the output of the \code{fitAll} function, which fits the GPA model for all possible pairs of GWAS datasets, as input.
}
\value{
}
\references{
Kortemeier E, Ramos PS, Hunt KJ, Kim HJ, Hardiman G, and Chung D (2017),
"ShinyGPA: An interactive and dynamic visualization toolkit for genetic studies."
}
\author{ Dongjun Chung, Emma Kortemeier }
\seealso{
\code{\link{fitAll}}.
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
