
.ff_emStep <- function( gwasPval, annMat, useAnn, pleiotropyH0, empiricalNull,
		maxIter, stopping, epsStopLL,
		binaryMat, betaAlpha, betaAlphaNull, pis, q1,
		lbPi, lbBetaAlpha, lbQ, vDigit, verbose ) {
		
	out <- .Call( "cppEM",
		gwasPval, annMat, useAnn, pleiotropyH0, empiricalNull,
		maxIter, stopping, epsStopLL,
		binaryMat, betaAlpha, betaAlphaNull, pis, q1,
		lbPi, lbBetaAlpha, lbQ, vDigit, verbose, PACKAGE = "GPA" )
	return(out)
}
