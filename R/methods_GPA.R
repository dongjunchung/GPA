
# generic methods for "GPA" class

setMethod(
    f="show",
    signature="GPA",
    definition=function( object ) {    
        # summary of GPA fit
		
		vDigit <- 1000
		
        # constants
		
		nBin <- nrow(object@gwasPval)
		nGWAS <- ncol(object@gwasPval)
		nAnn <- ncol(object@annMat)
		
		binaryList <- vector( "list", nGWAS )
		for ( k in 1:nGWAS ) {
			binaryList[[k]] <- c( 0, 1 )
		}
		binaryMat <- expand.grid( binaryList )
		
		nComp <- nrow(binaryMat)
		combVec <- apply( binaryMat, 1, function(bm) paste( bm, collapse="" ) )
		
        # parameters for GPA
		
		emSetting <- object@setting
		
		pis <- object@fit$pis
		betaAlpha <- object@fit$betaAlpha
		if ( emSetting$empiricalNull ) {
			betaAlphaNull <- object@fit$betaAlphaNull
		}
		q1 <- object@fit$q1
		
		nPis <- length(pis)
		nAlpha <- length(betaAlpha)
		
		# estimate covariance matrix
		
		covMat <- cov( object, silent=TRUE )
		seVec <- sqrt(diag(covMat))
		
		# SE for pi00, using Delta method
		
		piSE <- seVec[ 1:(nPis-1) ]			
		gderiv <- as.matrix( c( rep( -1, (nPis-1) ), rep( 0, nrow(covMat) - (nPis-1) ) ) )
		pi00SE <- sqrt( t(gderiv) %*% covMat %*% gderiv )
		
		piSE <- c( pi00SE, piSE )
		
		betaAlphaSE <- seVec[ (nPis-1+1):(nPis-1+nAlpha) ]		
		if ( emSetting$empiricalNull ) {
			betaAlphaNullSE <- seVec[ (nPis-1+nAlpha+1):(nPis-1+nAlpha+nAlpha) ]		
		}
		
		# output
        
        cat( "Summary: GPA model fitting results (class: GPA)\n" )
        cat( "--------------------------------------------------\n" )
        cat( "Data summary:\n" )
        cat( "\tNumber of GWAS data: ", nGWAS , "\n", sep="" )
		cat( "\tNumber of SNPs: ", nBin , "\n", sep="" )
        if ( emSetting$useAnn ) {       
			cat( "\tNumber of annotation data: ", nAnn, "\n", sep="" )
		} else {
			cat( "\tNumber of annotation data: (not provided)\n", sep="" )
		}
		cat( "Model setting:\n" )
        if ( emSetting$empiricalNull ) {       
			cat( "\tNull distribution is empirically estimated.\n" )
		} else {
			cat( "\tTheoretical null distribution is assumed.\n" )
		}
        if ( emSetting$pleiotropyH0 ) {       
			cat( "\tGPA is fitted under H0 of pleiotropy LRT.\n" )
		}
		cat( "Parameter estimates (standard errors):\n" )
		cat( "\talpha: ", paste( round(betaAlpha*vDigit)/vDigit, collapse=" " ), "\n", sep="" )
		cat( "\t     ( ", paste( round(betaAlphaSE*vDigit)/vDigit, collapse=" " ), " )", "\n", sep="" )
        if ( emSetting$empiricalNull ) {  
			cat( "\talphaNull: ", paste( round(betaAlphaNull*vDigit)/vDigit, collapse=" " ), "\n", sep="" )
			cat( "\t  ( ", paste( round(betaAlphaNullSE*vDigit)/vDigit, collapse=" " ), " )", "\n", sep="" )
		}
		cat( "\tGWAS combination: ", paste( combVec, collapse=" " ), "\n", sep="" )
		cat( "\tpi: ", paste( round(pis*vDigit)/vDigit, collapse=" " ), "\n", sep="" )
		cat( "\t  ( ", paste( round(piSE*vDigit)/vDigit, collapse=" " ), " )", "\n", sep="" )
		if ( emSetting$useAnn ) {
			cat( "\tq:\n" )
			for ( d in 1:nAnn ) { 
				cat( "\tAnnotation #",d,":\n", sep="" )			
				cat( "\t    ", paste( round(q1[d,]*vDigit)/vDigit, collapse=" " ), "\n", sep="" )
				
				if ( emSetting$empiricalNull ) {
					locQ1 <- (nPis-1+2*nAlpha+nPis*(d-1)+1):(nPis-1+2*nAlpha+nPis*d)
				} else {
					locQ1 <- (nPis-1+nAlpha+nPis*(d-1)+1):(nPis-1+nAlpha+nPis*d)
				}
				cat( "\t  ( ", paste( round(seVec[locQ1]*vDigit)/vDigit, collapse=" " ), " )", "\n", sep="" )
			}
		}
		if ( emSetting$useAnn ) {
			q1ratio <- q1ratioSE <- matrix( NA, nrow(q1), (ncol(q1)-1) )
			for ( d in 1:nAnn ) {
				# estimates
				
				q1ratio[d,] <- q1[d,-1] / q1[d,1]
				
				# SE
				
				for ( j in 2:ncol(q1) ) {
					qderiv <- rep( 0, nrow(q1)*ncol(q1) )
					qderiv[ ncol(q1) * (d-1) + 1 ] <- - q1[d,j] / q1[d,1]^2
					qderiv[ ncol(q1) * (d-1) + j ] <- 1 / q1[d,1]
					
					gderiv <- as.matrix( c( rep( 0, nrow(covMat) - nrow(q1)*ncol(q1) ), qderiv ) )
					q1ratioSE[d,(j-1)] <- sqrt( t(gderiv) %*% covMat %*% gderiv )
				}			
			}
			
			cat( "\n" )
			cat( "\tRatio of q over baseline (",combVec[1],"):\n", sep="" )
			cat( "\tGWAS combination: ", paste( combVec[-1], collapse=" " ),"\n", sep="" )
			for ( d in 1:nAnn ) {    	
				cat( "\tAnnotation #",d,":\n", sep="" )
				cat( "\t    ", paste( round(q1ratio[d,]*vDigit)/vDigit, collapse=" " ), "\n", sep="" )
				cat( "\t  ( ", paste( round(q1ratioSE[d,]*vDigit)/vDigit, collapse=" " ), " )", "\n", sep="" )
			}
		}
        cat( "--------------------------------------------------\n" )
    }
)

setMethod(
    f="print",
    signature="GPA",
    definition=function( x ) {
        # return posterior probability matrix
		
		return(x@fit$Z)
    }
)

setMethod(
    f="fdr",
    signature="GPA",
    definition=function( object, pattern=NULL ) {
      if ( is.null(pattern) ) {
        # if pattern is not specified, return marginal local FDR
		
		margfdr <- 1 - object@fit$Zmarg
		colnames(margfdr) <- colnames(object@gwasPval)
		
		return(margfdr)
      } else {
        # if pattern is provided, return FDR matching the pattern
        
        ppmat <- print(object)
        
        # parse pattern & convert it to regular expression
        
        ptvec <- strsplit( pattern, "" )[[1]]  
        if ( length(grep( "(1|\\*)", ptvec )) < length(ptvec) ) {
          stop( "Invalid 'pattern' argument! Please use only '1' or '*' for the pattern." )
        }
        ptreg <- paste( gsub( "\\*", "(0|1)", ptvec ), collapse="" )
        locp <- grep( ptreg, colnames(ppmat) )
        
        # return related FDR
        
        if ( length(locp) > 0 ) {
          psel <- 1 - rowSums(ppmat[ , locp, drop=FALSE ])
          return(psel)
        } else {
          stop( "Invalid 'pattern' argument! Please provide correct pattern." )
        }
      }
    }
)

setMethod(
    f="cov",
    signature="GPA",
    definition=function( object, silent=FALSE, vDigitEst=1000, vDigitSE=1000, ... ) {
        # calculate covariance matrix using empirical information matrix
		
		.covGPA( object, silent=silent, vDigitEst=vDigitEst, vDigitSE=vDigitSE )
    }
)

setMethod(
    f="se",
    signature="GPA",
    definition=function( object, ... ) {
        # return estimates & standard error
		
		# covariance matrix
		
		covEst <- .covGPA( object, silent=TRUE, vDigitEst=1000, vDigitSE=1000 )			
		seVec <- list()
	
		# extract estimates
		
		empiricalNull <- object@setting$empiricalNull
		useAnn <- object@setting$useAnn
		
		pis <- object@fit$pis
		betaAlpha <- object@fit$betaAlpha
		if ( empiricalNull ) {
			betaAlphaNull <- object@fit$betaAlphaNull
		}
		
		nGWAS <- length(betaAlpha)
		
		if ( useAnn ) {
			q1 <- object@fit$q1
			nAnn <- nrow(q1)
		}
		
		# extract information from estimates
					
		binaryList <- vector( "list", nGWAS )
		for ( k in 1:nGWAS ) {
			binaryList[[k]] <- c( 0, 1 )
		}
		binaryMat <- expand.grid( binaryList )
		
		nComp <- nrow(binaryMat)
		combVec <- apply( binaryMat, 1, function(bm) paste( bm, collapse="" ) )
		
		nPis <- length(pis)
		nAlpha <- length(betaAlpha)
		
		if ( useAnn ) {
			nQ1 <- nAnn * nComp
		}
		
		# SE for pi, using Delta method
		
		piSE <- sqrt(diag(covEst))[ 1:(nPis-1) ]			
		gderiv <- as.matrix( c( rep( -1, (nPis-1) ), 
			rep( 0, nrow(covEst) - (nPis-1) ) ) )
		pi00SE <- sqrt( t(gderiv) %*% covEst %*% gderiv )
		
		piSE <- c( pi00SE, piSE )
		names(piSE)[1] <- paste("pi_",names(pis)[1],sep="")
		
		# return estimates and SE	
		
		locPi <- 1:(nPis-1)
		locAlpha <- (nPis-1+1):(nPis-1+nAlpha)		
		
		seVec$betaAlpha <- sqrt(diag(covEst))[locAlpha]
		
		if ( empiricalNull ) {
			locAlpha0 <- (nPis-1+nAlpha+1):(nPis-1+nAlpha+nAlpha)			
			seVec$betaAlphaNull <- sqrt(diag(covEst))[locAlpha0]
		}
		
		seVec$pis <- piSE
		
		if ( useAnn ) {    
			# q
			
			seVec$q1 <- c()			
			for ( d in 1:nAnn ) {
				if ( empiricalNull ) {
					locQ1 <- (nPis-1+2*nAlpha+nPis*(d-1)+1):(nPis-1+2*nAlpha+nPis*d)
				} else {
					locQ1 <- (nPis-1+nAlpha+nPis*(d-1)+1):(nPis-1+nAlpha+nPis*d)
				}
				
				seVec$q1 <- rbind( seVec$q1, sqrt(diag(covEst))[locQ1] )
			}
			
			# ratio of q
			
			q1ratioSE <- matrix( NA, nrow(q1), (ncol(q1)-1) )
			for ( d in 1:nAnn ) {				
				for ( j in 2:ncol(q1) ) {
					qderiv <- rep( 0, nrow(q1)*ncol(q1) )
					qderiv[ ncol(q1) * (d-1) + 1 ] <- - q1[d,j] / q1[d,1]^2
					qderiv[ ncol(q1) * (d-1) + j ] <- 1 / q1[d,1]
					
					gderiv <- as.matrix( c( rep( 0, nrow(covEst) - nrow(q1)*ncol(q1) ), qderiv ) )
					q1ratioSE[d,(j-1)] <- sqrt( t(gderiv) %*% covEst %*% gderiv )
				}			
			}			
			
			seVec$q1ratio <- q1ratioSE
		}		
		
		return( seVec )
    }
)

setMethod(
    f="estimates",
    signature="GPA",
    definition=function( object, ... ) {
        # return parameter estimates
		
		paramEst <- list()
		
		paramEst$pis <- object@fit$pis
		paramEst$betaAlpha <- object@fit$betaAlpha
		if ( object@setting$empiricalNull ) { 
			paramEst$betaAlphaNull <- object@fit$betaAlphaNull
		}
		if ( object@setting$useAnn ) { 
			paramEst$q1 <- object@fit$q1
			
			# ratio of q
			
			paramEst$q1ratio <- matrix( NA, nrow(paramEst$q1), (ncol(paramEst$q1)-1) )
			for ( d in 1:nrow(paramEst$q1) ) {			
				paramEst$q1ratio[d,] <- paramEst$q1[d,-1] / paramEst$q1[d,1]
			}	
		}
		
		return(paramEst)
    }
)
