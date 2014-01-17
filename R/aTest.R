
aTest <- function( fitWithoutAnn, fitWithAnn, vDigit=1000 ) {
	
	empiricalNull <- fitWithoutAnn@setting$empiricalNull
	
	pis <- fitWithoutAnn@fit$pis
	betaAlpha <- fitWithoutAnn@fit$betaAlpha
	if ( empiricalNull ) {
		betaAlphaNull <- fitWithoutAnn@fit$betaAlphaNull
	}
	
	# constant
	
	nBin <- nrow(fitWithAnn@gwasPval)	
	nGWAS <- ncol(fitWithAnn@gwasPval)	
	
	binaryList <- vector( "list", nGWAS )
	for ( k in 1:nGWAS ) {
		binaryList[[k]] <- c( 0, 1 )
	}
	binaryMat <- expand.grid( binaryList )
	
	nComp <- nrow(binaryMat)
	combVec <- apply( binaryMat, 1, function(bm) paste( bm, collapse="" ) )
	
	#annMat <- as.matrix(fitWithAnn@annMat)	
	nAnn <- ncol(fitWithAnn@annMat)
		
	nPis <- length(pis)
	nAlpha <- length(betaAlpha)
	
	# LRT: numerator
	
	q1 <- matrix( colSums(fitWithAnn@annMat) / nrow(fitWithAnn@annMat), nAnn, nComp )
		
	betaDist <- matrix( NA, nBin, nComp )		
	for ( k in 1:nGWAS ) {
		betaDist[,k] <- dbeta( fitWithAnn@gwasPval[,k], betaAlpha[k], 1 )
	}
	if ( empiricalNull ) {
		betaNullDist <- matrix( NA, nBin, nComp )		
		
		for ( k in 1:nGWAS ) {
			betaNullDist[,k] <- dbeta( fitWithAnn@gwasPval[,k], betaAlphaNull[k], 1 )
		}
	}
    
    llmat <- matrix( NA, nBin, nComp )	
        
    for ( g in 1:nComp ) {
        # mixing proportion
	    
	    llmat[,g] <- pis[g]
        
        # emission for GWAS	   
	    
	    for ( k in 1:nGWAS ) {
	        if( binaryMat[g,k] == 1 ) {
	        	# signal SNP		        	
	        	     
	        	llmat[,g] <- llmat[,g] * betaDist[,k]
        	} else if ( empiricalNull ) {
	        	# null SNP, if null is empirically estimated
	        	
	        	llmat[,g] <- llmat[,g] * betaNullDist[,k]
        	}
	    }
	        	
    	# emission for annotation data
	
    	for ( d in 1:nAnn ) {
        	llmat[,g] <- llmat[,g] * 
        		c( ( 1 - q1[d,g] ), q1[d,g] )[ ( fitWithAnn@annMat[,d] + 1 ) ]
    	}
    }
    
	ll0 <- sum( log( rowSums( llmat ) ) )
	
	# LRT: denominator
	
	ll1 <- fitWithAnn@fit$loglik[ length(fitWithAnn@fit$loglik) ]
	
	# LRT
	
	LRT <- -2 * ( ll0 - ll1 )
	
	# calculate p-value
	
	degree.of.freedom <- nAnn * ( nComp - 1 )	
	pvalue <- pchisq( LRT, degree.of.freedom, lower.tail=FALSE )
	
	# SE for q1
	
	cov.GPA.wAnn <- cov( fitWithAnn, silent=TRUE )
	seVec <- sqrt(diag(cov.GPA.wAnn))
	
	# summary	
    
	message( "Hypothesis testing for annotation enrichment" )
	message( "( Note: This version of test is designed for single annotation data )" )
	message( "--------------------------------------------------" )
	message( "GWAS combination: ", paste( combVec, collapse=" " ) )
	message( "q1:" )
	for ( d in 1:nAnn ) {    	
    	message( "    ", paste( round(fitWithAnn@fit$q1[d,]*vDigit)/vDigit, collapse=" " ) )
    	
    	if ( empiricalNull ) {
    		locQ1 <- (nPis-1+2*nAlpha+nPis*(d-1)+1):(nPis-1+2*nAlpha+nPis*d)
		} else {
			locQ1 <- (nPis-1+nAlpha+nPis*(d-1)+1):(nPis-1+nAlpha+nPis*d)
		}
    	message( "  ( ", paste( round(seVec[locQ1]*vDigit)/vDigit, collapse=" " ), " )" )
	}
	message( " " )
    message( "test statistics: ", paste( round(LRT*vDigit)/vDigit, collapse=" " ) )
    message( "p-value: ", pvalue, collapse=" " )
	message( "--------------------------------------------------" )
	
	return( list( q=fitWithAnn@fit$q1, statistics=LRT, pvalue=pvalue ) )
}
