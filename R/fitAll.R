
# fit GPA for all possible combination

fitAll <- function( pmat, 
  maxIter=2000, stopping="relative", epsStopLL=1e-10, 
  parallel=FALSE, nCore=8 ) {

  # consider all possible cases
  
  combs <- t(combn( 1:ncol(pmat), 2 ))
  comblist <- split( combs, 1:nrow(combs) )
  
  
  # fit GPA
  
  message("-------------------------------------------------------------")
  message( "Fitting GPA for all possible combination of GWAS datasets...", sep="" )
  message("-------------------------------------------------------------")
  
  if ( parallel==TRUE ) {
    fit.list.GPA <- mclapply( comblist, function(cc) {
      message( "\tGWAS pair: ",cc[1]," and ",cc[2], sep="" )
      fit.GPA.cc <- GPA( pmat[ , c(cc[1],cc[2]) ], NULL, pleiotropyH0=FALSE, 
        maxIter=maxIter, stopping=stopping, epsStopLL=epsStopLL, verbose=0 )
      return(fit.GPA.cc)
    }, mc.cores=nCore )
  } else {
    fit.list.GPA <- lapply( comblist, function(cc) {
      message( "\tGWAS pair: ",cc[1]," and ",cc[2], sep="" )
      fit.GPA.cc <- GPA( pmat[ , c(cc[1],cc[2]) ], NULL, pleiotropyH0=FALSE, 
        maxIter=maxIter, stopping=stopping, epsStopLL=epsStopLL, verbose=0 )
      return(fit.GPA.cc)
    })
  }
  
  message("")
  message("Done!")
  message("")
  
  
  # fit GPA under H0
  
  message("-------------------------------------------------------------")
  message( "Fitting GPA for all possible combination under H0...", sep="" )
  message("-------------------------------------------------------------")
  
  if ( parallel==TRUE ) {
    fit.list.GPA.H0 <- mclapply( comblist, function(cc) {
      message( "Fit GPA for the GWAS ",cc[1]," and ",cc[2], sep="" )
      fit.GPA.H0.cc <- GPA( pmat[ , c(cc[1],cc[2]) ], NULL, pleiotropyH0=TRUE, 
        maxIter=maxIter, stopping=stopping, epsStopLL=epsStopLL, verbose=0 )
      return(fit.GPA.H0.cc)
    }, mc.cores=nCore )
  } else {
    fit.list.GPA.H0 <- lapply( comblist, function(cc) {
      message( "Fit GPA for the GWAS ",cc[1]," and ",cc[2], sep="" )
      fit.GPA.H0.cc <- GPA( pmat[ , c(cc[1],cc[2]) ], NULL, pleiotropyH0=TRUE, 
        maxIter=maxIter, stopping=stopping, epsStopLL=epsStopLL, verbose=0 )
      return(fit.GPA.H0.cc)
    })
  }
  
  message("")
  message("Done!")
  message("")
  
  
  # calculate pleiotropy test p-values
  
  message("-------------------------------------------------------------")
  message( "Calculating pleiotropy test p-values...", sep="" )
  message("-------------------------------------------------------------")

  plist <- lapply( 1:length(fit.list.GPA), function(i) {
  	ll.H1 <- fit.list.GPA[[i]]@fit$loglik[ length(fit.list.GPA[[i]]@fit$loglik) ]
  	ll.H0 <- fit.list.GPA.H0[[i]]@fit$loglik[ length(fit.list.GPA.H0[[i]]@fit$loglik) ]
  	
  	LRT <- -2 * ( ll.H0 - ll.H1 )
  	pval <- pchisq( LRT, 1, lower.tail=FALSE )
    return(pval)
  } )
  
  pvec <- unlist(plist)
  pTestPval <- matrix( 0, ncol=ncol(pmat), nrow=ncol(pmat) )
  rownames(pTestPval) <- colnames(pmat)
  colnames(pTestPval) <- colnames(pmat)
  
  for ( i in 1:nrow(combs) ) {
    pTestPval[ combs[i,1], combs[i,2] ] <- pvec[i]
    pTestPval[ combs[i,2], combs[i,1] ] <- pvec[i]
  }
  
  message("Done!")
  
  return(list( pmat=pmat, combs=combs, combList=comblist, 
    pTestPval=pTestPval, fitGPA=fit.list.GPA, fitH0=fit.list.GPA.H0 ))
}
