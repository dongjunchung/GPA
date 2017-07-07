
.covGPA <- function( object, silent=FALSE, vDigitEst=1000, vDigitSE=1000 ) {
	# check correctness of arguments
	
	if ( silent != TRUE & silent != FALSE ) {
		stop( "Inappropriate value for 'silent' argument. It should be either TRUE or FALSE." )
	}
	
	if ( vDigitEst %% 10 != 0 | vDigitEst <= 0 ) {
		stop( "Inappropriate value for 'vDigitEst' argument. It should be multiples of 10, e.g., 10, 100, ..." )
	}
	
	if ( vDigitSE %% 10 != 0 | vDigitSE <= 0 ) {
		stop( "Inappropriate value for 'vDigitSE' argument. It should be multiples of 10, e.g., 10, 100, ..." )
	}
	
	# load fits
	
	emSetting <- object@setting
	gwasPval <- object@gwasPval
	
	# convert p-value & annotation vector to matrix, if needed, for consistency
	
	if ( is.vector(gwasPval) ) {
		gwasPval <- as.matrix(gwasPval)
	}
	if ( emSetting$useAnn ) {   
		annMat <- object@annMat
	} else {
		annMat <- NULL
	}
	
	# extract estimates
	
	empiricalNull <- emSetting$empiricalNull
	Z <- object@fit$Z
	pis <- object@fit$pis
	betaAlpha <- object@fit$betaAlpha
	if ( empiricalNull ) {
		betaAlphaNull <- object@fit$betaAlphaNull
	}
	
	nBin <- nrow(Z)
	nGWAS <- length(betaAlpha)
	
	if ( !is.null(annMat) ) {
		q1 <- object@fit$q1
		nAnn <- ncol(annMat)
	}
	
	# extract information from estimates
		
	sumZ <- colSums(Z)	
	
	binaryList <- vector( "list", nGWAS )
	for ( k in 1:nGWAS ) {
		binaryList[[k]] <- c( 0, 1 )
	}
	binaryMat <- expand.grid( binaryList )
	
	nComp <- nrow(binaryMat)
	combVec <- apply( binaryMat, 1, function(bm) paste( bm, collapse="" ) )
	
	nPis <- length(pis)
	nAlpha <- length(betaAlpha)
	
	if ( !is.null(annMat) ) {
		nQ1 <- nAnn * nComp
	}
	
	# indexing
	
	locZero <- which( rowSums(binaryMat) == 0 )
	
	locNonzero <- 1:nPis
	locNonzero <- locNonzero[ -locZero ]
	
	namePis <- names(pis)[ -locZero ]
	
	########################################################################
    #                                                                      #
    #                 empirical observed information                       #
    #                                                                      #
    ########################################################################
	
    if ( empiricalNull ) {

	    if ( !is.null(annMat) ) {
	    	smat <- matrix( NA, nBin, nPis - 1 + 2 * nAlpha + nQ1 )
		} else {
			smat <- matrix( NA, nBin, nPis - 1 + 2 * nAlpha )
		}
		
	} else {
		
	    if ( !is.null(annMat) ) {
	    	smat <- matrix( NA, nBin, nPis - 1 + nAlpha + nQ1 )
		} else {
			smat <- matrix( NA, nBin, nPis - 1 + nAlpha )
		}
		
	}
	
	# pi
	
	for ( i in 1:(nPis-1) ) {
		smat[,i] <- Z[ , locNonzero[i] ] / pis[ locNonzero[i] ] -
			Z[ , locZero ] / pis[ locZero ]
	}
	
	# alpha
	
	for ( i in 1:nAlpha ) {
		smat[ , nPis-1+i ] <- ( log(gwasPval[,i]) + 1 / betaAlpha[i] ) * 
			rowSums( Z[ , binaryMat[,i] == 1, drop=FALSE ] )
	}
	
	if ( empiricalNull ) {
		
		# alpha0
		
		for ( i in 1:nAlpha ) {
			smat[ , nPis-1+nAlpha+i ] <- ( log(gwasPval[,i]) + 1 / betaAlphaNull[i] ) * 
				rowSums( Z[ , binaryMat[,i] == 0, drop=FALSE ] )
		}
		
		# q1
		
		if ( !is.null(annMat) ) {	
			k <- 1
			
			for ( d in 1:nAnn ) {
				for ( j in 1:nComp ) {
					smat[ , nPis-1+2*nAlpha+k ] <-
						Z[,j] * ( annMat[,d] / q1[d] - ( 1 - annMat[,d] ) / ( 1 - q1[d] ) )
					
					k <- k + 1
				}
			}
		}
		
	} else {
		
		# q1
		
		if ( !is.null(annMat) ) {	
			k <- 1
			
			for ( d in 1:nAnn ) {
				for ( j in 1:nComp ) {
					smat[ , nPis-1+nAlpha+k ] <-
						Z[,j] * ( annMat[,d] / q1[d] - ( 1 - annMat[,d] ) / ( 1 - q1[d] ) )
					
					k <- k + 1
				}
			}
		}
		
	}
		
	# empirical observed information
	
	infoEmp <- t(smat) %*% smat
	covEst <- solve(infoEmp)
	
	if ( empiricalNull ) {
	
		if ( !is.null(annMat) ) {	
			rownames(infoEmp) <- colnames(infoEmp) <- rownames(covEst) <- colnames(covEst) <-
				c( paste("pi_",namePis,sep=""), 
				paste("alpha_",1:nAlpha,sep=""), paste("alpha0_",1:nAlpha,sep=""),
				paste("q1_",rep(1:nAnn,each=nComp),"_",rep(1:nComp,nAnn),sep="") )
		} else {
			rownames(infoEmp) <- colnames(infoEmp) <- rownames(covEst) <- colnames(covEst) <-
				c( paste("pi_",namePis,sep=""), 
				paste("alpha_",1:nAlpha,sep=""), paste("alpha0_",1:nAlpha,sep="") )
		}
		
	} else {
		
		if ( !is.null(annMat) ) {	
			rownames(infoEmp) <- colnames(infoEmp) <- rownames(covEst) <- colnames(covEst) <-
				c( paste("pi_",namePis,sep=""), paste("alpha_",1:nAlpha,sep=""),
				paste("q1_",rep(1:nAnn,each=nComp),"_",rep(1:nComp,nAnn),sep="") )
		} else {
			rownames(infoEmp) <- colnames(infoEmp) <- rownames(covEst) <- colnames(covEst) <-
				c( paste("pi_",namePis,sep=""), paste("alpha_",1:nAlpha,sep="") )
		}
	
	}
		
	# SE for pi, using Delta method
	
	piSE <- sqrt(diag(covEst))[ 1:(nPis-1) ]			
	gderiv <- as.matrix( c( rep( -1, (nPis-1) ), 
		rep( 0, nrow(covEst) - (nPis-1) ) ) )
	pi00SE <- sqrt( t(gderiv) %*% covEst %*% gderiv )
	
	piSE <- c( pi00SE, piSE )
	
	# report
	
	if( !silent ) {
		locPi <- 1:(nPis-1)
		locAlpha <- (nPis-1+1):(nPis-1+nAlpha)		
		    
		message( " " )
		message( "alpha: ", paste( round(betaAlpha*vDigitEst)/vDigitEst, collapse=" " ) )
		message( "     ( ", paste( round(sqrt(diag(covEst))[locAlpha]*vDigitSE)/vDigitSE, collapse=" " ), " )" )
		
		if ( empiricalNull ) {
			locAlpha0 <- (nPis-1+nAlpha+1):(nPis-1+nAlpha+nAlpha)
			message( "alpha0: ", paste( round(betaAlphaNull*vDigitEst)/vDigitEst, collapse=" " ) )
			message( "     ( ", paste( round(sqrt(diag(covEst))[locAlpha0]*vDigitSE)/vDigitSE, collapse=" " ), " )" )
		}
		message( "GWAS combination: ", paste( combVec, collapse=" " ) )
		message( "pi: ", paste( round(pis*vDigitEst)/vDigitEst, collapse=" " ) )
		message( "  ( ", paste( round(piSE*vDigitSE)/vDigitSE, collapse=" " ), " )" )
		
		if ( !is.null(annMat) ) {    
			# q
			
			message( "q: " )
			for ( d in 1:nAnn ) {
				message( "Annotation #",d,":")
				if ( empiricalNull ) {
					locQ1 <- (nPis-1+2*nAlpha+nPis*(d-1)+1):(nPis-1+2*nAlpha+nPis*d)
				} else {
					locQ1 <- (nPis-1+nAlpha+nPis*(d-1)+1):(nPis-1+nAlpha+nPis*d)
				}
				
				message( "\t    ", paste( round(q1[d,]*vDigitEst)/vDigitEst, collapse=" " ) )
				message( "\t  ( ", paste( round(sqrt(diag(covEst))[locQ1]*vDigitSE)/vDigitSE, collapse=" " ), " )" )
			}
			
			# ratio of q
			
			q1ratio <- q1ratioSE <- matrix( NA, nrow(q1), (ncol(q1)-1) )
			for ( d in 1:nAnn ) {
				# estimates
				
				q1ratio[d,] <- q1[d,-1] / q1[d,1]
				
				# SE
				
				for ( j in 2:ncol(q1) ) {
					qderiv <- rep( 0, nrow(q1)*ncol(q1) )
					qderiv[ ncol(q1) * (d-1) + 1 ] <- - q1[d,j] / q1[d,1]^2
					qderiv[ ncol(q1) * (d-1) + j ] <- 1 / q1[d,1]
					
					gderiv <- as.matrix( c( rep( 0, nrow(covEst) - nrow(q1)*ncol(q1) ), qderiv ) )
					q1ratioSE[d,(j-1)] <- sqrt( t(gderiv) %*% covEst %*% gderiv )
				}			
			}
			
			message( " " )
			message( "Ratio of q over baseline (",combVec[1],"):" )
			message( "GWAS combination: ", paste( combVec[-1], collapse=" " ) )
			for ( d in 1:nAnn ) {
				message( "Annotation #",d,":")
				message( "\t    ", paste( round(q1ratio[d,]*vDigitEst)/vDigitEst, collapse=" " ) )
				message( "\t  ( ", paste( round(q1ratioSE[d,]*vDigitSE)/vDigitSE, collapse=" " ), " )" )
			}
		}
	}
	
	return( covEst )	
}