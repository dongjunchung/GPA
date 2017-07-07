
# main GPA function

GPA <- function( gwasPval, annMat=NULL, pleiotropyH0=FALSE, empiricalNull=FALSE,
	maxIter=2000, stopping="relative", epsStopLL=1e-10,
	initBetaAlpha=0.1, initPi=0.1, initQ=0.75,
	lbPi=NA, lbBetaAlpha=0.001, lbQ=0.001, lbPval=1e-30,
	vDigit=1000, verbose=1 ) {

	# check correctness of arguments
	
	if ( pleiotropyH0 != TRUE & pleiotropyH0 != FALSE ) {
		stop( "Inappropriate value for 'pleiotropyH0' argument. It should be either TRUE or FALSE." )
	}
	
	if ( empiricalNull != TRUE & empiricalNull != FALSE ) {
		stop( "Inappropriate value for 'empiricalNull' argument. It should be either TRUE or FALSE." )
	}
	
	if ( maxIter %% 1 != 0 | maxIter <= 0 ) {
		stop( "Inappropriate value for 'maxIter' argument. It should be positive integer." )
	}
	
	if ( stopping != "absolute" & stopping != "relative" & stopping != "aitken" ) {
		stop( "Inappropriate value for 'stopping' argument. It should be one of 'absolute' (absolute difference in log likelihood), 'relative' (relative difference in log likelihood), and 'aitken' (Aitken acceleration-based stopping rule)." )
	}
	
	if ( epsStopLL < 0 ) {
		stop( "Inappropriate value for 'epsStopLL' argument. It should be nonnegative." )
	}
	
	if ( initBetaAlpha <= 0 | initBetaAlpha >= 1 ) {
		stop( "Inappropriate value for 'initBetaAlpha' argument. It should be between zero and one." )
	}
	
	if ( initPi <= 0 | initPi >= 1 ) {
		stop( "Inappropriate value for 'initPi' argument. It should be between zero and one." )
	}
	
	if ( initQ <= 0 | initQ >= 1 ) {
		stop( "Inappropriate value for 'initQ' argument. It should be between zero and one." )
	}
	
	if ( !is.na(lbPi) && ( lbPi <= 0 | lbPi >= 1 ) ) {
		stop( "Inappropriate value for 'lbPi' argument. It should be NA ('lbPi' = 1 / [number of SNPs]) or some value between zero and one." )
	}
	
	if ( lbBetaAlpha <= 0 | lbBetaAlpha >= 1 ) {
		stop( "Inappropriate value for 'lbBetaAlpha' argument. It should be between zero and one." )
	}
	
	if ( lbQ <= 0 | lbQ >= 1 ) {
		stop( "Inappropriate value for 'lbQ' argument. It should be between zero and one." )
	}
	
	if ( lbPval <= 0 | lbPval >= 1 ) {
		stop( "Inappropriate value for 'lbPval' argument. It should be between zero and one." )
	}
	
	if ( vDigit %% 10 != 0 | vDigit <= 0 ) {
		stop( "Inappropriate value for 'vDigit' argument. It should be multiples of 10, e.g., 10, 100, ..." )
	}	
	
	if ( verbose != 0 & verbose != 1 & verbose != 2 & verbose != 3 ) {
		stop( "Inappropriate value for 'verbose' argument. It should be one of 0, 1, 2, or 3." )
	}
	
	# convert p-value & annotation vector to matrix, if needed, for consistency
	
	if ( !is.matrix(gwasPval) ) {
		gwasPval <- as.matrix(gwasPval)
	}
	if ( !is.null(annMat) ) {
		if ( !is.matrix(annMat) ) {
			annMat <- as.matrix(annMat)
		}
	}	
	
	# check correctness of data
	# gwasPval: M * K matrix (p-value, [ 0, 1 ] )
	# annMat: M * D matrix (binary, { 0, 1 } )
	
	if ( !is.null(annMat) ) {
		if ( nrow(gwasPval) != nrow(annMat) ) {
			stop( "Number of SNPs are different between p-value and annotation matrices. They should coincide. Please check your p-value and annotation matrices.")
		}
	}
	
	if ( any( gwasPval < 0 | gwasPval > 1 ) ) {
		stop( "Some p-values are smaller than zero or larger than one. p-value should be ranged between zero and one. Please check your p-value matrix." )
	}
	
	if ( !is.null(annMat) ) {
		if ( any( annMat != 0 & annMat != 1 ) ) {
			stop( "Some elements in annotation matrix has values other than zero or one. All the elements in annotation matrix should be either zero (not annotated) or one (annotated). Please check your annotation matrix." )
		}
	}
	
	# initialization of constants
	
	nBin <- nrow(gwasPval)
	nGWAS <- ncol(gwasPval)
	
	if ( !is.null(annMat) ) {
		nAnn <- ncol(annMat)
	}	
	
	# report setting
	
	message( "Info: Number of GWAS data: ", nGWAS )
	if ( !is.null(annMat) ) {
		message( "Info: Number of annotation data: ", nAnn )
	}
	
	if ( empiricalNull ) {
		message( "Info: Null distribution is empirically estimated." )
	} else {
		message( "Info: Theoretically null distribution is used." )
	}
	
	if ( !is.null(annMat) ) {
		message( "Info: Annotation data is provided." )
		message( "Info: SNPs will be prioritized using annotation data." )
	} else {
		message( "Info: No annotation data is provided." )
	}
	
	if ( pleiotropyH0 ) {
		message( "Info: Fit the GPA model under H0 of the pleitropy test." )
	}
	
	if ( is.na(lbPi) ) {
		message( "Info: Lower bound for pi estimates is set to 1 / [number of SNPs]." )
		lbPi <- 1 / nBin
	}

	# set zero p-values to small values to avoid log(0)
	
	if ( verbose >= 1 ) {
		if ( any( gwasPval < lbPval ) ) {
			message( "Info: Some SNPs have p-values close to zero." )
			message( "Info: Number of SNPs with p-values close to zero: ", length(which( gwasPval < lbPval )) )
			message( "Info: p-values for these SNPs are set to ", lbPval )
			
			gwasPval[ gwasPval < lbPval ] <- lbPval
		}
	}
	
	# define binary matrix for multiple GWAS datas
	
	binaryList <- vector( "list", nGWAS )
	for ( k in 1:nGWAS ) {
		binaryList[[k]] <- c( 0, 1 )
	}
	binaryMat <- as.matrix(expand.grid( binaryList ))
	
	nComp <- nrow(binaryMat)
	combVec <- apply( binaryMat, 1, function(bm) paste( bm, collapse="" ) )
	
	# initialization of parameters
	
	#betaAlpha <- rep( initBetaMean / ( 1 - initBetaMean ), nGWAS )
	betaAlpha <- rep( initBetaAlpha, nGWAS )
	betaAlphaNull <- rep( 1, nGWAS )
	
	pis <- initPi^rowSums(binaryMat)
	if ( sum(pis[2:nComp]) >= 1 ) {
		pis[2:nComp] <- 0.4 * pis[2:nComp] / sum(pis[2:nComp])
	}
	pis[ rowSums(binaryMat) == 0 ] <- 1 - sum( pis[ rowSums(binaryMat) != 0 ] )
	if ( any( pis < lbPi ) ) {
		pis[ pis < lbPi ] <- lbPi
	}
	pis <- pis / sum(pis)
	
	if ( !is.null(annMat) ) {
		q1 <- matrix( initQ, nAnn, nComp )
	} else {
		q1 <- matrix( rep(1,4), 2, 2 )
	}
    
    # EM step
    
    useAnn <- as.numeric( !is.null(annMat) )
    if ( is.null(annMat) ) {
	    annMat <- matrix( rep(1,4), 2, 2 )
    }
    stoppingCode <- switch( stopping,
    	absolute = 1,
    	relative = 2,
    	aitken = 3
    )
	    
    emResult <- .ff_emStep( 
    	gwasPval = gwasPval, annMat = annMat, useAnn = useAnn, 
		pleiotropyH0 = as.numeric( pleiotropyH0 ), empiricalNull = as.numeric( empiricalNull ),
		maxIter = maxIter, stopping = stoppingCode, epsStopLL = epsStopLL, 		
		binaryMat = binaryMat, betaAlpha = betaAlpha, betaAlphaNull = betaAlphaNull, pis = pis, q1 = q1, 
    	lbPi = lbPi, lbBetaAlpha = lbBetaAlpha, lbQ = lbQ, vDigit = vDigit, verbose = verbose )
	
	# post processing of EM fitting
	
	#emResult$empiricalNull <- empiricalNull
	
	emResult$finalIter <- emResult$maxIter
			
	emResult$piMat <- emResult$piMat[ 1:(emResult$finalIter), , drop=FALSE ]
	emResult$betaAlphaMat <- emResult$betaAlphaMat[ 1:(emResult$finalIter), , drop=FALSE ]
	emResult$betaAlphaNullMat <- emResult$betaAlphaNullMat[ 1:(emResult$finalIter), , drop=FALSE ]
	emResult$loglik <- emResult$loglik[ 1:(emResult$finalIter) ]
	emResult$q1Array <- emResult$q1Array[ 1:(emResult$finalIter), , drop=FALSE ]
	
	colnames(emResult$Z) <- names(emResult$pis) <- colnames(emResult$piMat) <- combVec
	
	rownames(emResult$piMat) <- rownames(emResult$betaAlphaMat) <- names(emResult$loglik) <-
		paste( "iteration_", 1:nrow(emResult$betaAlphaMat), sep="" )
	if ( empiricalNull ) {
		rownames(emResult$betaAlphaNullMat) <- 
			paste( "iteration_", 1:nrow(emResult$betaAlphaMat), sep="" )
	}
		
	if ( useAnn ) {
		rownames(emResult$q1) <- colnames(annMat)			
		colnames(emResult$q1) <- combVec
	} else {
		emResult$q1 <- NA
		emResult$q1Mat <- matrix( NA, 1, 1 )
	}
	
	# summarizing setting for GPA
	
	emSetting <- list()

	emSetting$useAnn <- useAnn
	emSetting$pleiotropyH0 <- pleiotropyH0
	emSetting$empiricalNull <- empiricalNull
		
	emSetting$maxIter <- maxIter
	emSetting$stopping <- stopping
	emSetting$epsStopLL <- epsStopLL
	emSetting$initBetaAlpha <- initBetaAlpha
	emSetting$initPi <- initPi
	emSetting$initQ1 <- initQ
	emSetting$lbPi <- lbPi
	emSetting$lbBetaAlpha <- lbBetaAlpha
	emSetting$lbQ <- lbQ
	emSetting$lbPval <- lbPval
	
	# return object by creating GPA class object
	
    new( "GPA",
        fit = emResult, setting = emSetting, gwasPval = gwasPval, annMat = annMat )
}
