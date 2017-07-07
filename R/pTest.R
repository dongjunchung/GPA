
pTest <- function( fit, fitH0, vDigit=1000 ) {
	# check correctness of arguments
	
	if ( !is( fit, "GPA" ) ) {
		stop( " Input for 'fit' argument is not 'GPA' class object. Please check the input." )
	}
	
	if ( !is( fitH0, "GPA" ) ) {
		stop( " Input for 'fitH0' argument is not 'GPA' class object. Please check the input." )
	}
	
	if ( vDigit %% 10 != 0 | vDigit <= 0 ) {
		stop( "Inappropriate value for 'vDigit' argument. It should be multiples of 10, e.g., 10, 100, ..." )
	}
	
	# load fits
	
	pis <- fit@fit$pis
	covMat <- cov( fit, silent=TRUE )
	
	# constant
	
	nGWAS <- length(fit@fit$betaAlpha)	
	
	binaryList <- vector( "list", nGWAS )
	for ( k in 1:nGWAS ) {
		binaryList[[k]] <- c( 0, 1 )
	}
	binaryMat <- expand.grid( binaryList )
	
	nComp <- nrow(binaryMat)
	combVec <- apply( binaryMat, 1, function(bm) paste( bm, collapse="" ) )
	
	# SE for pi, using Delta method
	
	piSE <- sqrt(diag(covMat))[ 1:(length(pis)-1) ]			
	gderiv <- as.matrix( c( rep( -1, (length(pis)-1) ), 
		rep( 0, nrow(covMat) - (length(pis)-1) ) ) )
	pi00SE <- sqrt( t(gderiv) %*% covMat %*% gderiv )
	
	piSE <- c( pi00SE, piSE )
	
	# calculate LRT & p-value
	
	ll.H1 <- fit@fit$loglik[ length(fit@fit$loglik) ]
	ll.H0 <- fitH0@fit$loglik[ length(fitH0@fit$loglik) ]
	
	LRT <- -2 * ( ll.H0 - ll.H1 )
	pvalue <- pchisq( LRT, 1, lower.tail=FALSE )
	
	# summary	
    
	cat( "Hypothesis testing for pleiotropy\n" )
	cat( "--------------------------------------------------\n" )
	cat( "GWAS combination: ", paste( combVec, collapse=" " ), "\n" )
    cat( "pi: ", paste( round(pis*vDigit)/vDigit, collapse=" " ), "\n" )
    cat( "  ( ", paste( round(piSE*vDigit)/vDigit, collapse=" " ), " )\n" )
	cat( "\n" )
    cat( "test statistics: ", paste( round(LRT*vDigit)/vDigit, collapse=" " ), "\n" )
    cat( "p-value: ", pvalue, "\n", collapse=" " )
	cat( "--------------------------------------------------\n" )
	
	return( list( pi=pis, piSE=piSE, statistics=LRT, pvalue=pvalue ) )
}
