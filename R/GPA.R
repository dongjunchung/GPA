
# main GPA function

GPA <- function( gwasPval, annMat=NULL, pleiotropyH0=FALSE, empiricalNull=FALSE,
	maxIter=2000, stopping="relative", epsStopLL=1e-10,
	initBetaMean=0.1, initPi=0.1, initQ1=0.75,
	lbPi=NA, lbBetaAlpha=0.001, lbQ=0.001, lbPval=1e-30,
	vDigit=1000, verbose=1 ) {
		
	# gwasPval: M * K matrix (p-value, [ 0, 1 ] )
	# annMat: M * D matrix (binary, { 0, 1 } )
	
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
	
	if ( !is.null(annMat) ) {
		if ( nrow(gwasPval) != nrow(annMat) ) {
			stop( "Number of SNPs are different between p-value and annotation matrices. Please check your p-value and annotation matrices!")
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
	
	betaAlpha <- rep( initBetaMean / ( 1 - initBetaMean ), nGWAS )
	betaAlphaNull <- rep( 1, nGWAS )
	
	pis <- initPi^rowSums(binaryMat)
	if ( any( pis < lbPi ) ) {
		pis[ pis < lbPi ] <- lbPi
	}
	pis[ rowSums(binaryMat) == 0 ] <- 1 - sum( pis[ rowSums(binaryMat) != 0 ] )
	
	if ( !is.null(annMat) ) {
		q1 <- matrix( initQ1, nAnn, nComp )
	} else {
		q1 <- as.matrix(1)
	}
    
    # EM step
    
    useAnn <- as.numeric( !is.null(annMat) )
    if ( is.null(annMat) ) {
	    annMat <- as.matrix(1)
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
	emResult$piMat = as.matrix( emResult$piMat[ 1:(emResult$finalIter), ] )
	emResult$betaAlphaMat = as.matrix( emResult$betaAlphaMat[ 1:(emResult$finalIter), ] )
	emResult$betaAlphaNullMat = as.matrix( emResult$betaAlphaNullMat[ 1:(emResult$finalIter), ] )
	emResult$loglik = emResult$loglik[ 1:(emResult$finalIter) ]
	emResult$q1Array = as.matrix( emResult$q1Array[ 1:(emResult$finalIter), ] )
	
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
	emSetting$initBetaMean <- initBetaMean
	emSetting$initPi <- initPi
	emSetting$initQ1 <- initQ1
	emSetting$lbPi <- lbPi
	emSetting$lbBetaAlpha <- lbBetaAlpha
	emSetting$lbQ <- lbQ
	emSetting$lbPval <- lbPval
	
	# return object by creating GPA class object
	
    new( "GPA",
        fit = emResult, setting = emSetting, gwasPval = gwasPval, annMat = annMat )
}
