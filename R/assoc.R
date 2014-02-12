
setMethod(
    f="assoc",
    signature="GPA",
    definition=function( object, FDR=0.05, fdrControl="global" ) {
		
		# check arguments
		
		if ( fdrControl != "global" & fdrControl != "local" ) {
			stop( "Invalid value for 'fdrControl' argument! It should be either 'global' or 'local'." )
		}
		
		if ( FDR < 0 | FDR > 1 ) {
			stop( "Invalid value for 'FDR' argument! It should be between zero and one." )
		}
		
		# association mapping
		
		fdrmat <- fdr(object)
		amat <- matrix( 0, nrow(fdrmat), ncol(fdrmat) )
		
		if ( fdrControl == "local" ) {
			# local FDR control
			
			message( "Info: Association mapping based on the local FDR control at level ", FDR, "." )
			
			amat[ fdrmat <= FDR ] <- 1
		} else if ( fdrControl == "global" ) {
			# global FDR control
			
			message( "Info: Association mapping based on the global FDR control at level ", FDR, "." )
		
			# direct approach for FDR control
			
			for ( j in 1:ncol(amat) ) {
				pp <- fdrmat[,j]
				pp.ordered <- sort(pp)
				pp.cum <- cumsum( pp.ordered ) / c(1:length(pp))
				cutoff <- max( pp.ordered[ pp.cum <= FDR ] )
				amat[ pp <= cutoff, j ] <- 1
			}
			
		}
		
		return(amat)
	}
)
